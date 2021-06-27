
% 
% VLFeat  License:
% Copyright (C) 2007-11, Andrea Vedaldi and Brian Fulkerson
% Copyright (C) 2012-13, The VLFeat Team
% All rights reserved.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function [pos2] = find_matches(I1,pos1,I2)
% Installing external open-source ibrary install
 run('vlfeat-0.9.21\toolbox\vl_setup.m')
load stereoPointPairs

global width

%For some reason, if the images were kept as double precision, i was
%getting errors
I1 = im2uint8(I1);
I2 = im2uint8(I2);
%convert images to grayscale double 
i1 = single(rgb2gray(I1));
width = size(i1, 2);

i2 = single(rgb2gray(I2));

% Threshold for choosing points of interest
siftThreshold = 2;   %feeling smart, might delete later

% Finding points of interest
[fa, da] = vl_sift(i1, 'PeakThresh', siftThreshold) ;
[fb, db] = vl_sift(i2, 'PeakThresh', siftThreshold) ;

% Creating list of putative matches
[matches, scores] = vl_ubcmatch(da, db) ;

fb(1,:) = fb(1,:) + size(i1,2) ;

% Coordinates of keypoints
points1 = [fa(1,:);fa(2,:)];
points2 = [fb(1,:);fb(2,:)];


% Create dataset with corresponding keypoints according to SIFT matches
%FIrst column represents x1, second y1, third x2, fourth y2
[rows,cols] = size(matches);
data = zeros(cols,4);
for i = 1:cols
    data(i, 1:2) = points1(:,matches(1,i));
    data(i, 3:4) = points2(:, matches(2, i));
end

%Running the Ransac algorithm with the Matlab built in funtion
%By using custom implementations of the functions that fit the model and
%evaluate distaces
[model, inliers] = ransac(data, @fitTransformationMatrix, @evaluateDistanceFromModel, 20, 10, 'MaxNumTrials', 2000);

%After the model is evaluated, we use it to find the corresponding
%locations to the entries of pos1(similar to the way positions are
%calculated when evaluating the distance)
matched = [];
[rows,cols] = size(pos1);
for i = 1 : rows
    initLoc = [pos1(i, 1); pos1(i, 2); 1];
    unscaled = model * initLoc;
    resLoc = [unscaled(1)/ unscaled(3), unscaled(2)/ unscaled(3)];
    matched = [matched; resLoc];
end

%converting to actual coordinates, not coordinates when images are put next
%to each other
matched(:,1) = matched(:, 1) - width;

%print matches
%used for testing
% figure;
% showMatchedFeatures(I1,I2,pos1,matched,'montage','PlotOptions',{'ro','go','b--'});
% title('Point matches');

pos2 = matched;

end


function [tMatrix] = fitTransformationMatrix(data)
%Computing Homography matrix
matrix = [];
[rows,cols] = size(data);
for i = 1:rows
    matrix = [matrix; 
                   -1 * data(i, 1),  -1 *data(i, 2), -1, 0, 0, 0, data(i, 3) * data(i, 1),  data(i, 3) * data(i, 4),data(i, 3);
                    0, 0, 0, -1 *  data(i, 1), -1 *data(i, 2), -1,data(i, 4) * data(i, 1) data(i, 4) * data(i, 2), data(i, 4)];
end

[U, S, V] = svd(matrix);
h = V(:,9);
tMatrix = reshape(h,3,3)';
end

function [distances] = evaluateDistanceFromModel(model, data)
%Using the Homography matrix to calcukate new coordinates and calculate
%their distance from the actual match
global width
distances = [];
[rows, cols] = size(data);
for i = 1:rows
    origin = [data(i, 1); data(i, 2); 1];
    actualMatch = [data(i, 3), data(i, 4)];
    temp = model*origin;
    t1 = temp(1)/temp(3);
    t2 = temp(2)/temp(3);
    if t1 < width
        curr = 9999;
    else 
        curr = norm(actualMatch - [abs(t1), abs(t2)]);
    end
    
    distances = [distances,curr];
end
distances;
end






