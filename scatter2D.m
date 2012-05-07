try
	aa = load("AA")
	scatter(aa(:,1), aa(:,2), 'r')
catch
	disp ("no aa cluster")
end_try_catch

hold

try
	ab = load("AB")
	scatter(ab(:,1), ab(:,2), 'g')
catch
	disp ("no ab cluster")
end_try_catch

try 
	bb = load("BB")
	scatter(bb(:,1), bb(:,2), 'b')
catch
	disp ("no bb cluster")
end_try_catch

hold('off')
