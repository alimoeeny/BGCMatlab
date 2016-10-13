function response = Run_one_rds(itterations, stim_disparity, filtersize, lc, ls, rc, rs)
%
%	response = Run_one_rds(itterations, stim_disparity, rf_disparity, filtersize, sd, freq)
%
%	returns 4 by itterations array of responses of 4 rf's to random contrast stimuli
%	[left_eye_cos, left_eye_sin, right_eye_cos, right_eye_sin]
%	... where the rf's are odd and even gabors shifted the appropriate disparity
%
%
% Setup_odd_gabor RFs returns sin,cos for left, -sin,cos for right, for
% odd symetric tuning
	response = [0 0 0 0];
% calculate response to different random 1d patterns
	for i = 1:itterations
		[leftstim, rightstim] = Stim_rds_1d(stim_disparity, filtersize);
		a(1) = sum(leftstim .* lc);
		a(2) = sum(leftstim .* ls);
		a(3) = sum(rightstim .* rc);
		a(4) = sum(rightstim .* rs);
		if(i == 1)
		  response = a;
		else
		  response = [response;a];
		end
	end
