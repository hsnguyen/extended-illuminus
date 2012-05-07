try
	aa = load("AA")
	scatter3(aa(:,1), aa(:,2), aa(:,3),'r')
catch
	disp ("no aa cluster")
end_try_catch

hold ('on')

try
	ab = load("AB")
	scatter3(ab(:,1), ab(:,2), ab(:,3), 'g')
catch
	disp ("no ab cluster")
end_try_catch

try 
	bb = load("BB")
	scatter3(bb(:,1), bb(:,2), bb(:,3), 'b')
catch
	disp ("no bb cluster")
end_try_catch

hold('off')
