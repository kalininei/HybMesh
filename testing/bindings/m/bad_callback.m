function ret=bad_callback(s1, s2, p1, p2)
	out = sprintf('%s %f --- %s %f\n', s1, p1, s2, p2);
	disp(out);
	if p1 < 0.7
		ret = 0;
	else
		ret = 1;
	end
end
