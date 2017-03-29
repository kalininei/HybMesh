function ret=callback(s1, s2, p1, p2)
	out = sprintf('%s %f --- %s %f\n', s1, p1, s2, p2);
	disp(out);
	ret = 0;
end
