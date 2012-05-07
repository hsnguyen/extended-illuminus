try
	aa = load("AA")
#	for i=1:length(aa(:,1))
#		[aa(i,1), aa(i,2)] = detrans(aa(i,1), aa(i,2))
#	end
	scatter(aa(:,1), aa(:,2), 'r')
catch
	disp ("no aa cluster")
end_try_catch

hold

try
	ab = load("AB")
#	for i=1:length(ab(:,1))
#                [ab(i,1), ab(i,2)] = detrans(ab(i,1), ab(i,2))
#        end
	scatter(ab(:,1), ab(:,2), 'g')
catch
	disp ("no ab cluster")
end_try_catch

try 
	bb = load("BB")
#	for i=1:length(bb(:,1))
#                [bb(i,1), bb(i,2)] = detrans(bb(i,1), bb(i,2))
#        end
	scatter(bb(:,1), bb(:,2), 'b')
catch
	disp ("no bb cluster")
end_try_catch

length1 = length(bb(:,1))
length2 = length(ab(:,1))
length3 = length(aa(:,1))

scatter(aa(length3, 1), aa(length3, 2), 'k', '+')
scatter(ab(length2, 1), ab(length2, 2), 'k', '+')
scatter(bb(length1, 1), bb(length1, 2), 'k', '+')

hold('off')
