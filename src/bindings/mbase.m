classdef Worker
	properties
		connection;
		callback = @(n1, n2, p1, p2) 1;
	end
	methods(Access={?HybMesh})
		%low level communication
		function ret=require_connection(self, server_path)
			%TODO
		end
		function ret=get_signal(self)
			%TODO
		end
		function [sz, data] = get_data(self)
			%TODO
		end
		function send_signal(self, sig)
			%TODO
		end
		function send_data(self, sz, data)
			%TODO
		end
		function break_connection(self)
			%TODO
		end
	end
	methods (Access=private)
		function _send_command(self, func, com){
			self.send_data(size(func, 2), func);
			self.send_data(size(com, 2), com);
			self.send_signal("C");
		end
		function ret=_wait_for_signal(self)
			ret=self.get_signal(self)
		end
		function ret= _read_buffer(self)
			[ret, ~] = self.get_data();
		end
		function _apply_callback(self, buf)
			p = typecast(buf(1:16), "double");
			lens = typecast(buf(16:24), "int32");
			s1 = buf(24:24+lens(1));
			s2 = buf(24+lens(1):24+lens(1)+lens(2));
			r = self.callback(s1, s2, p(1), p(2));
			if r == 0
				self.send_signal("G");
			else
				self.send_signal("S");
			end
		end
		function ret=__parse_vecstring(self, str)
		end
	end
	methods
		function ret=_apply_command(self, func, com)
			self._send_command(input);
			while (1)
				sig = self._wait_for_signal();
				switch (sig)
					%normal return
					case 'R'
						ret = self._read_buffer();
						return;
					%callback
					case 'B'
						self._apply_callback(self._read_buffer());
					%exception return
					case 'I'
						error("User interrupted");
					case 'E'
						error(self._read_buffer());
					%something went wrong
					case '0'
						error("Server stopped working");
					otherwise
						error("Invalid client instruction");
				end
			end
		end
		% string -> cpp type
		function ret=_to_int(self, str)
			ret = str2num(str);
		end
		function ret=_to_double(self, str)
			ret = str2num(str);
		end
		function ret=_to_point(self, str)
			if (strcmp(str, "None"))
				ret = nan;
			else
				ret = str2num(str);
			end
		end
		function ret=_to_grid(self, str)
			str = self.__parse_vecstring(str)
			ret = Grid2D(str(1), self);
		end
		function ret=_to_cont(self, str)
			str = self.__parse_vecstring(str)
			ret = Contour2D(str(1), self);
		end
		function ret=_to_grid3(self, str)
			str = self.__parse_vecstring(str)
			ret = Grid3D(str(1), self);
		end
		function ret=_to_vecint(self, str)
			ret = str2num(str)
		end
		function ret=_to_veccont(self, str)
			str = self.__parse_vecstring(str)
			ret(size(str, 1)) = Contour2D;
			for i = [1:size(str, 1)]
				ret(i) = Contour2D(str(i), self);
			end
		end
		function ret=_to_vecgrid(self, str)
			str = self.__parse_vecstring(str)
			ret(size(str, 1)) = Grid2D;
			for i = [1:size(str, 1)]
				ret(i) = Grid2D(str(i), self);
			end
		end
		function ret=_to_vecsurface(self, str)
			str = self.__parse_vecstring(str)
			ret(size(str, 1)) = Surface3D;
			for i = [1:size(str, 1)]
				ret(i) = Surface3D(str(i), self);
			end
		end
		function ret=_to_vecgrid3(self, str)
			str = self.__parse_vecstring(str)
			ret(size(str, 1)) = Grid3D;
			for i = [1:size(str, 1)]
				ret(i) = Grid3D(str(i), self);
			end
		end
		% cpp type -> string
		function ret=_tos_bool(self, val)
			if val
				ret = "True";
			else
				ret = "False";
			end
		end
		function ret=_tos_int(self, val)
			ret =  sprintf('%i', val);
		end
		function ret=_tos_double(self, val)
			ret = sprintf('%.16g', val);
		end
		function ret=_tos_point(self, val)
			if (isnan(val))
				ret = "None";
			else
				ret = sprintf('[%.16g, %.16g]', val(1), val(2));
			end
		end
		function ret=_tos_point3(self, val)
			if (isnan(val))
				ret = "None";
			else
				ret = sprintf('[%.16g, %.16g, %.16g]', val(1), val(2), val(3));
			end
		end
		function ret=_tos_vecbyte(self, val)
			ret = val;
		end
		function ret=_tos_vecint(self, val)
			as = sprintf('%i, ', vec(:));
			ret = strcat('[', as, ']');
		end
		function ret=_tos_vecdouble(self, val)
			as = sprintf('%.16g, ', vec(:));
			ret = strcat('[', as, ']');
		end
		function ret=_tos_vecstring(self, val)
			as = sprintf("\"%s\", ", val(:));
			ret = strcat('[', as, ']')
		end
		function ret=_tos_vecpoint(self, val)
			ret = sprintf("[%f, %f], ", val(1, :), val(2, :))
			ret = strcat('[', ret, ']');
		end
		function ret=_tos_object(self, val);
			if (isnan(val) || size(val.sid, 2) == 0)
				ret = "None"
			else
				ret = '"' + val.sid + '"';
			end
		end
		function ret=_tos_vecobject(self, val);
			ret = sprintf("\"%s\", ", val(:).sid);
			ret = strcat('[', ret, ']');
		end
		% vecbyte -> cpp typ
		function ret=_to_map_int_double_raw(self, val)
			sz = typecast(val(1:4), "int32");
			ret = zeros(2, sz);
			pos = 5;
			for i=[1:sz]
				ret(1, i) = typecast(val(pos:pos+4), "int32"); 
				pos = pos + 4;
				ret(2, i) = typecast(val(pos:pos+8), "double"); 
				pos = pos + 8;
			end
		end
		function ret=_to_vecdouble_raw(self, val)
			sz = typecast(val(1:4), "int32");
			ret = zeros(1, sz); 
			pos = 5;
			for i=[1:sz]
				ret(1, i) = typecast(val(pos:pos+8), "double");
				pos = pos + 8;
			end
		end
	end
end

classdef Contour2D
	methods(Access={?HybMesh})
		worker;
		sid;
	end
	methods
		function self=Contour2D(sid="", worker=0)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Contour2D
	end
end

classdef Grid2D
	methods(Access={?HybMesh})
		worker;
		sid;
	end
	methods
		function self=Grid2D(sid="", worker=0)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Grid2D
	end
end

classdef Surface3D
	methods(Access={?HybMesh})
		worker;
		sid;
	end
	methods
		function self=Surface3D(sid="", worker=0)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Surface3D
	end
end

classdef Grid3D
	methods(Access={?HybMesh})
		worker;
		sid;
	end
	methods
		function self=Grid3D(sid="", worker=0)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Grid3D
	end
end

classdef Hybmesh < handle
	properties (Access=private)
		Worker worker;
	end

	methods
		function ret=Hybmesh(path)
			ret.worker = Worker;
			ret.worker.connection = ret.worker.require_connection(path);
		end
		function assign_callback(self, cb)
			ret.worker.callback = cb;
		end
		function delete(self)
			if self.worker.connection != -1
				self.worker.break_connection();
				self.worker.connection = -1;
			end
		end

		%>>$Hybmesh
	end
end
