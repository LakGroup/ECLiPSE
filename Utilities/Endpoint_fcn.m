function is_end_point = Endpoint_fcn(nhood)
% -------------------------------------------------------------------------
% A function that makes a lookup table that is able to determine whether or
% not a pixel in an image is an endpoint or not.
% Example on how to use it:
%   EndPointLUT = makelut(@Endpoint_fcn, 3);
% -------------------------------------------------------------------------
% Input:
%   nhood:  Size of the LUT neighbourhood (2 or 3 in general).
%
% Output:
%   is_end_point:   Logical value for whether a pixel is an endpoint or
%                   not.
% -------------------------------------------------------------------------
% Code obtained from:
%   Digital Image Processing using MATLAB (ISBN: 978-0-9820854-0-0) by R.C.
%   Gonzalez, R.E. Woods & S.L. Eddins (code can be found on page 508).
% Contact:
%   siewert.hugelier@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   xxx
% -------------------------------------------------------------------------

% Define the hit-or-miss interval matrices for end points (obtained from
% Gonzalez 7 Woods, 2008; ISBN: 978-0-13-168728-8)
interval1 = [0 1 0; -1 1 -1; -1 -1 -1];
interval2 = [1 -1 -1; -1 1 -1; -1 -1 -1];

% The end points are defined by the hit-or-miss interval matrices or any of
% their 90 degrees rotations.
for i = 1:4
    C = bwhitmiss(nhood, rot90(interval1,i)); % Interval 1 + rotations.
    D = bwhitmiss(nhood, rot90(interval2,i)); % Interval 2 + rotations.

    % Check if the neighbourhood matches any of the endpoint
    % configurations. If yes, then return a true valie and finish the
    % function.
    if (C(2,2) == 1) || (D(2,2) == 1)
        is_end_point = true;
        return
    end
end

% If the pixel neighbourhood matched none of the endpoint configurations,
% then return a false value.
is_end_point = false;