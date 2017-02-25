classdef Worker
	methods
		function ret=_apply_command(self, func, com)
		end
		% string -> cpp type
		function ret=_to_int(self, str)
		end
		function ret=_to_double(self, str)
		end
		function ret=_to_point(self, str)
		end
		function ret=_to_grid(self, str)
		end
		function ret=_to_cont(self, str)
		end
		function ret=_to_grid3(self, str)
		end
		function ret=_to_vecint(self, str)
		end
		function ret=_to_veccont(self, str)
		end
		function ret=_to_vecgrid(self, str)
		end
		function ret=_to_vecsurface(self, str)
		end
		function ret=_to_vecgrid3(self, str)
		end
		% cpp type -> string
		function ret=_tos_bool(self, val)
		end
		function ret=_tos_int(self, val)
		end
		function ret=_tos_double(self, val)
		end
		function ret=_tos_point(self, val)
		end
		function ret=_tos_point3(self, val)
		end
		function ret=_tos_vecbyte(self, val)
		end
		function ret=_tos_vecint(self, val)
		end
		function ret=_tos_vecdouble(self, val)
		end
		function ret=_tos_vecstring(self, val)
		end
		function ret=_tos_vecpoint(self, val)
		end
		function ret=_tos_vecobject(self, val);
		end
		% vecbyte -> cpp typ
		function ret=_to_map_int_double_raw(self, val)
		end
		function ret=_to_vecdouble_raw(self, val)
		end
	end
end

classdef Contour2D
	properties
		worker;
		sid;
	end
	methods
		function self=Contour2D(sid, worker)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Contour2D
	end
end

classdef Grid2D
	properties
		worker;
		sid;
	end
	methods
		function self=Grid2D(sid, worker)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Grid2D
	end
end

classdef Surface3D
	properties
		worker;
		sid;
	end
	methods
		function self=Surface3D(sid, worker)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Surface3D
	end
end

classdef Grid3D
	properties
		worker;
		sid;
	end
	methods
		function self=Grid3D(sid, worker)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Grid3D
	end
end

classdef Hybmesh
	properties
		Worker worker;
	end
	methods
		%>>$Hybmesh
	end
end
