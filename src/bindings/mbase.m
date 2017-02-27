classdef HybmeshWorker
	properties
		connection;
		callback = @(n1, n2, p1, p2) 1;
	end
	methods(Access={?HybMesh})
		%low level communication
		function ret=require_connection(self, server_path)
			ret = core_hmconnection_oct(1, server_path);
		end
		function ret=get_signal(self)
			ret = core_hmconnection_oct(2, self.connection);
		end
		function data = get_data(self)
			ret = core_hmconnection_oct(4, self.connection);
		end
		function send_signal(self, sig)
			core_hmconnection_oct(3, self.connection, sig);
		end
		function send_data(self, data)
			core_hmconnection_oct(5, self.connection, data);
		end
		function break_connection(self)
			core_hmconnection_oct(6, self.connection);
		end
	end
	methods (Access=private)
		function _send_command(self, func, com){
			self.send_data(func);
			self.send_data(com);
			self.send_signal("C");
		end
		function ret=_wait_for_signal(self)
			ret=self.get_signal(self)
		end
		function ret= _read_buffer(self)
			ret = self.get_data();
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
			ret = HybmeshGrid2D(str(1), self);
		end
		function ret=_to_cont(self, str)
			str = self.__parse_vecstring(str)
			ret = HybmeshContour2D(str(1), self);
		end
		function ret=_to_grid3(self, str)
			str = self.__parse_vecstring(str)
			ret = HybmeshGrid3D(str(1), self);
		end
		function ret=_to_vecint(self, str)
			ret = str2num(str)
		end
		function ret=_to_veccont(self, str)
			str = self.__parse_vecstring(str)
			ret(size(str, 1)) = HybmeshContour2D;
			for i = [1:size(str, 1)]
				ret(i) = HybmeshContour2D(str(i), self);
			end
		end
		function ret=_to_vecgrid(self, str)
			str = self.__parse_vecstring(str)
			ret(size(str, 1)) = HybmeshGrid2D;
			for i = [1:size(str, 1)]
				ret(i) = HybmeshGrid2D(str(i), self);
			end
		end
		function ret=_to_vecsurface(self, str)
			str = self.__parse_vecstring(str)
			ret(size(str, 1)) = HybmeshSurface3D;
			for i = [1:size(str, 1)]
				ret(i) = HybmeshSurface3D(str(i), self);
			end
		end
		function ret=_to_vecgrid3(self, str)
			str = self.__parse_vecstring(str)
			ret(size(str, 1)) = HybmeshGrid3D;
			for i = [1:size(str, 1)]
				ret(i) = HybmeshGrid3D(str(i), self);
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

classdef HybmeshContour2D
	methods(Access={?HybMesh})
		worker;
		sid;
	end
	methods
		function self=HybmeshContour2D(sid="", worker=0)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Contour2D
	end
end

classdef HybmeshGrid2D
	methods(Access={?HybMesh})
		worker;
		sid;
	end
	methods
		function self=HybmeshGrid2D(sid="", worker=0)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Grid2D
	end
end

classdef HybmeshSurface3D
	methods(Access={?HybMesh})
		worker;
		sid;
	end
	methods
		function self=HybmeshSurface3D(sid="", worker=0)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Surface3D
	end
end

classdef HybmeshGrid3D
	methods(Access={?HybMesh})
		worker;
		sid;
	end
	methods
		function self=HybmeshGrid3D(sid="", worker=0)
			self.sid = sid;
			self.worker = worker;
		end
		%>>$Grid3D
	end
end

classdef Hybmesh < handle
	properties (Access=private)
		HybmeshWorker worker;
	end
	methods(Static)
		function ret=hybmesh_exec_path(newpath)
			persistent path = "hybmesh.exe";  %>>$EXEPATH
			if nargin
				path = newpath;
			end
			ret = path;
		end
	end
	methods
		function ret=Hybmesh(path)
			ret.worker = HybmeshWorker;
			ret.worker.connection = ret.worker.require_connection(
				ret.hybmesh_exec_path());
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
