function [x, y] = detrans(c, s)
	x = exp(s) * (1 + c) / 2;
	y = exp(s) * (1 - c) / 2;
end
