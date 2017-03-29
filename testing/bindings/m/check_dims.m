function check_dims(v1, v2)
	if (length(v1) ~= length(v2)), error('check failed'); end
	for i = 1:length(v2)
		if (v1(i) ~= v2(i)), error('check failed'); end
	end
end


