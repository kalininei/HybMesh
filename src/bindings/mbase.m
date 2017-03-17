classdef HybmeshWorker < handle
	properties(Access={?Hybmesh})
		connection;
		callback = @(~, ~, ~, ~) 0;
	end
	methods(Access=private)
		%low level communication
		function ret=require_connection(~, server_path)
			ret = core_hmconnection_oct(1, server_path);
		end
		function ret=get_signal(self)
			ret = core_hmconnection_oct(2, self.connection);
		end
		function ret = get_data(self)
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
		function _send_command(self, func, com)
			self.send_signal('C');
			self.send_data(func);
			self.send_data(com);
		end
		function ret=_wait_for_signal(self)
			ret=self.get_signal(self);
		end
		function ret= _read_buffer(self)
			ret = self.get_data();
		end
		function _apply_callback(self, buf)
			p = typecast(buf(1:16), 'double');
			lens = typecast(buf(17:24), 'int32');
			s1 = buf(25:24+lens(1));
			s2 = buf(25+lens(1):24+lens(1)+lens(2));
			r = self.callback(s1, s2, p(1), p(2));
			if r == 0
				self.send_signal('G');
			else
				self.send_signal('S');
			end
		end

		% removes ending []; splits by ,; strips substrings from ''', ''', ' ';
		function ret=__parse_vecstring(~, str)
			s = strtrim(str);
			if (s(1)=='[' && s(end)==']')
				s = substr(s, 2, length(s)-2);
			end
			if (strcmp(s, '[]') || strcmp(s, 'None') || length(s) == 0)
				ret=[];
				return;
			end
			ret1 = {};
			for cit = strsplit(s, ',')
				it = char(cit);
				pos1 = 1;
				pos2 = length(it);
				while (pos1 <= length(it) &&
				       (it(pos1) == ' ' ||
				        it(pos1) == '''' ||
				        it(pos1) == '''')) pos1 = pos1 + 1;
				end
				while (pos2 >= pos1 &&
				       (it(pos2) == ' ' ||
				        it(pos2) == '''' ||
				        it(pos2) == '''')) pos2 = pos2 - 1;
				end
				if (pos2 >= pos1)
					ret1{end+1}  = substr(it, pos1, pos2-pos1+1);
				end
			end
			bracket_level = 0;
			ret = {};
			for cit=ret1
				it = char(cit);
				if (bracket_level == 0)
					ret{end+1} = it;
				else
					ret{end} = strcat(ret{end}, ',', it);
				end
				if (it(1) == '[') bracket_level = bracket_level+1; end
				if (it(end) == ']') bracket_level = bracket_level-1; end
			end
		end
		function ret = isnone(~, obj)
			if isnumeric(obj) && length(obj) == 1
				ret=isnan(obj);
			else
				ret=false;
			end
		end
	end
	methods(Access={?Hybmesh, ?HybmeshObject})
		function ret = HybmeshWorker(exepath)
			ret.connection=ret.require_connection(exepath);
		end
		function free(self)
			if self.connection ~= -1
				self.break_connection();
				self.connection = -1;
			end
		end
		function ret=_apply_command(self, func, com)
			self._send_command(func, com);
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
						error('Interrupted by user');
					case 'E'
						error(self._tos_vecbyte(self._read_buffer()));
					%something went wrong
					case '0'
						error('Server stopped working');
					otherwise
						error('Invalid client instruction');
				end
			end
		end
		function ret = _to_object(self, ret, str)
			ss = self.__parse_vecstring(str);
			if (length(ss) == 0)
				ret=nan;
			else
				if strcmp(ss{1}, 'None')
					ret = nan;
				else
					ret.sid = ss{1};
					ret.worker = self;
				end
			end
		end
		function ret = _to_vecobject(self, cls, str)
			ss = self.__parse_vecstring(str);
			if (length(ss) == 0)
				ret={};
				return;
			end
			for i=[1:length(ss)]
				ret{i} = self._to_object(cls, ss{i});
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
			if (strcmp(str, 'None'))
				ret = nan;
			else
				ret = str2num(str);
			end
		end
		function ret=_to_grid(self, str)
			ret = self._to_object(HybmeshGrid2D, str);
		end
		function ret=_to_cont(self, str)
			ret = self._to_object(HybmeshContour2D, str);
		end
		function ret=_to_grid3(self, str)
			ret = self._to_object(HybmeshGrid3D, str);
		end
		function ret=_to_surface(self, str)
			ret = self._to_object(HybmeshSurface3D, str);
		end
		function ret=_to_veccont(self, str)
			ret = self._to_vecobject(HybmeshContour2D, str);
		end
		function ret=_to_vecgrid(self, str)
			ret = self._to_vecobject(HybmeshGrid2D, str);
		end
		function ret=_to_vecsurface(self, str)
			ret = self._to_vecobject(HybmeshSurface3D, str);
		end
		function ret=_to_vecgrid3(self, str)
			ret = self._to_vecobject(HybmeshGrid3D, str);
		end
		function ret=_to_vecint(self, str)
			ret = str2num(str);
		end
		% cpp type -> string
		function ret=_tos_bool(self, val)
			if val
				ret = 'True';
			else
				ret = 'False';
			end
		end
		function ret=_tos_string(self, val)
			if self.isnone(val)
				ret = 'None'
			else
				ret = strcat('''', val, '''');
			end
		end
		function ret=_tos_int(self, val)
			ret =  sprintf('%i', val);
		end
		function ret=_tos_double(self, val)
			ret = sprintf('%.16g', val);
		end
		function ret=_tos_point(self, val)
			if (self.isnone(val))
				ret = 'None';
			else
				ret = sprintf('[%.16g, %.16g]', val);
			end
		end
		function ret=_tos_point3(self, val)
			if (self.isnone(val))
				ret = 'None';
			else
				ret = sprintf('[%.16g, %.16g, %.16g]', val);
			end
		end
		function ret=_tos_vecbyte(self, val)
			ret = val;
		end
		function ret=_tos_vecint(self, val)
			as = sprintf('%i, ', val);
			ret = strcat('[', as, ']');
		end
		function ret=_tos_vecdouble(self, val)
			as = sprintf('%.16g, ', val);
			ret = strcat('[', as, ']');
		end
		function ret=_tos_vecstring(self, val)
			if length(val) == 0
				ret = '[]';
			else
				as = sprintf('''%s'', ', val);
				ret = strcat('[', as, ']');
			end
		end
		function ret=_tos_vecpoint(self, val)
			if self.isnone(val)
				ret = 'None'
			elseif length(val) == 0
				ret = '[]';
			else
				ret = sprintf('[%.16g, %.16g], ', val');
				ret = strcat('[', ret, ']');
			end
		end
		function ret=_tos_object(self, val)
			if self.isnone(val)
				ret = 'None';
			else
				ret = self._tos_string(val.sid);
			end
		end
		function ret=_tos_vecobject(self, val)
			if self.isnone(val)
				ret = 'None';
			elseif length(val) == 0
				ret = '[]';
			else
				ret = strcat('''', val{1}.sid, '''');
				for i=2:length(val)
					ret = strcat(ret, ',''', val{i}.sid, '''');
				end
				ret = strcat('[', ret, ']');
			end
		end
		% vecbyte -> cpp typ
		function ret=_to_vec_int_double_raw(self, val)
			sz = typecast(val(1:4), 'int32');
			ret = zeros(sz, 2);
			pos = 5;
			for i=1:sz
				ret(i, 1) = typecast(val(pos:pos+3), 'int32'); 
				pos = pos + 4;
				ret(i, 2) = typecast(val(pos:pos+7), 'double'); 
				pos = pos + 8;
			end
		end
		function ret=_to_vecdouble_raw(self, val)
			ret = typecast(val(5:end), 'double');
		end
		function ret=_to_vecint_raw(self, val)
			ret = typecast(val(5:end), 'int32');
		end
	end
end

classdef HybmeshObject
	properties (Access={?HybmeshWorker, ?HybmeshObject})
		sid;
		worker;
	end
	methods
		%>>$ObjectA
	end
end

classdef HybmeshObject2D < HybmeshObject
	methods
		%>>$Object2D
	end
end

classdef HybmeshObject3D < HybmeshObject
	methods
		%>>$Object3D
	end
end

classdef HybmeshContour2D < HybmeshObject2D
	methods
		%>>$Contour2D
	end
end

classdef HybmeshGrid2D < HybmeshObject2D
	methods
		%>>$Grid2D
	end
end

classdef HybmeshSurface3D < HybmeshObject3D
	methods
		%>>$Surface3D
	end
end

classdef HybmeshGrid3D < HybmeshObject3D
	methods
		%>>$Grid3D
	end
end

classdef Hybmesh < handle
	properties (Access=private)
		worker;
	end
	methods(Static)
		function ret=hybmesh_exec_path(newpath)
			persistent path = 'hybmesh.exe';  %>>$EXEPATH
			if nargin
				path = newpath;
			end
			ret = path;
		end
		function ret=hybmesh_lib_path(newpath)
			persistent path = 'hybmesh.so';  %>>$LIBPATH
			if nargin
				path = newpath;
			end
			ret = path;
		end
	end
	methods
		function ret=Hybmesh()
			addpath(ret.hybmesh_lib_path());
			ret.worker = HybmeshWorker(ret.hybmesh_exec_path());
		end
		function assign_callback(self, cb)
		% Assigns callback function as @(string, string, double, doube)->int
		% which returns 1 for cancellation request and 0 otherwise.

			self.worker.callback = cb;
		end
		function reset_callback(self)
		% Sets default silent callback

			self.worker.callback = @(~,~,~,~) 0;
		end
		function delete(self)
			self.worker.free();
		end
		%>>$Hybmesh
	end
end
