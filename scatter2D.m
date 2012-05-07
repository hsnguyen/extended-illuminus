#function scatter2D(name)

#hold('off')
try
	aa = load("AA")
	length3 = length(aa(:,1))
#	for i=1:length(aa(:,1))
#		[aa(i,1), aa(i,2)] = detrans(aa(i,1), aa(i,2))
#	end
	scatter(aa(1:(length3-1),1), aa(1:(length3-1),2), 'r')
	hold('on')
	scatter(aa(length3, 1), aa(length3, 2), 'k', '+')
catch
	disp ("no aa cluster")
end_try_catch

hold('on')

try
	ab = load("AB")
	length2 = length(ab(:,1))
#	for i=1:length(ab(:,1))
#               [ab(i,1), ab(i,2)] = detrans(ab(i,1), ab(i,2))
#       end
	scatter(ab(1:(length2-1),1), ab(1:(length2-1),2), 'g')
	scatter(ab(length2, 1), ab(length2, 2), 'k', '+')
catch
	disp ("no ab cluster")
end_try_catch

try 
	bb = load("BB")
	length1 = length(bb(:,1))
#	for i=1:length(bb(:,1))
#                [bb(i,1), bb(i,2)] = detrans(bb(i,1), bb(i,2))
#        end
	scatter(bb(1:(length1-1),1), bb(1:(length1-1),2), 'b')
	scatter(bb(length1, 1), bb(length1, 2), 'k', '+')
catch
	disp ("no bb cluster")
end_try_catch

hold('off')
#end
