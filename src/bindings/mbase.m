classdef HybmeshCallerOctave
	properties (Access={?HybmeshWorker})
		connection;
	end
	methods (Access={?HybmeshWorker})
		function ret=HybmeshCallerOctave(exepath)
			ret.connection=ret.require_connection(exepath);
		end
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
end

classdef HybmeshCallerMatlab
	properties (Access={?HybmeshWorker})
		connection;
	end
	methods (Access={?HybmeshWorker})
		function ret=HybmeshCallerMatlab(exepath)
			if ~libisloaded('libcore_hmconnection_matlab')
				loadlibrary('libcore_hmconnection_matlab', ...
					'libcore_hmconnection_matlab.h');
			end
			ret.connection=ret.require_connection(exepath);
		end
		%low level communication
		function ret=require_connection(~, server_path)
			ret = calllib('libcore_hmconnection_matlab', ...
				'require_connection', server_path);
			ret = int32(ret);
		end
		function ret=get_signal(self)
			ret = calllib('libcore_hmconnection_matlab', ...
				'get_signal', self.connection);
			ret = cast(ret, 'char');
		end
		function ret = get_data(self)
			sz = calllib('libcore_hmconnection_matlab', ...
				'get_data1', self.connection);
			if sz == 0
				ret = [];
			else
				buf = libpointer('voidPtr', zeros(1, sz, 'uint8'));
				calllib('libcore_hmconnection_matlab', ...
					'get_data2', self.connection, int32(sz), buf);
				ret = buf.Value;
				clear buf;
			end
		end
		function send_signal(self, sig)
			calllib('libcore_hmconnection_matlab', ...
				'send_signal', self.connection, int8(sig));
		end
		function send_data(self, data)
			calllib('libcore_hmconnection_matlab', ...
				'send_data', self.connection, int32(length(data)), data);
		end
		function break_connection(self)
			calllib('libcore_hmconnection_matlab', ...
				'break_connection', self.connection);
		end
	end
end

classdef HybmeshWorker < handle
	properties(Access=private)
		caller;
	end
	properties(Access={?Hybmesh})
		callback = @(a, b, c, d) 0;
	end
	methods (Access=private)
		function send_command(self, func, com)
			self.caller.send_signal('C');
			self.caller.send_data(func);
			self.caller.send_data(com);
		end
		function ret=wait_for_signal(self)
			ret=self.caller.get_signal();
		end
		function ret= read_buffer(self)
			ret = self.caller.get_data();
		end
		function apply_callback(self, buf)
			p = typecast(buf(1:16), 'double');
			lens = typecast(buf(17:24), 'int32');
			s1 = buf(25:24+lens(1));
			s2 = buf(25+lens(1):24+lens(1)+lens(2));
			r = self.callback(s1, s2, p(1), p(2));
			if r == 0
				self.caller.send_signal('G');
			else
				self.caller.send_signal('S');
			end
		end

		% removes ending []; splits by ,; strips substrings from ''', ''', ' ';
		function ret=parse_vecstring(~, str)
			s = strtrim(str);
			if (s(1)=='[' && s(end)==']')
				s = s(2:end-1);
			end
			if (strcmp(s, '[]') || strcmp(s, 'None') || isempty(s))
				ret=[];
				return;
			end
			ret1 = {};
			strrange = strsplit(s, ',');
			ret1{length(strrange)} = [];
			iend = 0;
			for cit = strrange
				it = char(cit);
				pos1 = 1;
				pos2 = length(it);
				while (pos1 <= length(it) && ...
				       (it(pos1) == ' ' || ...
				        it(pos1) == '''' || ...
				        it(pos1) == '''')), pos1 = pos1 + 1;
				end
				while (pos2 >= pos1 && ...
				       (it(pos2) == ' ' || ...
				        it(pos2) == '''' || ...
				        it(pos2) == '''')), pos2 = pos2 - 1;
				end
				if (pos2 >= pos1)
					iend = iend + 1;
					ret1{iend} = it(pos1:pos2);
				end
			end
			ret = {};
			if iend == 0, return; end
			ret{iend} = [];
			bracket_level = 0;
			iend2 = 0;
			for cit=ret1
				it = char(cit);
				if (bracket_level == 0)
					iend2 = iend2 + 1;
					ret{iend2} = it;
				else
					ret{iend2} = strcat(ret{iend2}, ',', it);
				end
				if (it(1) == '['), bracket_level = bracket_level+1; end
				if (it(end) == ']'), bracket_level = bracket_level-1; end
			end
			ret(iend2+1:end) = [];
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
			try
				% try to load matlab procedure
				ret.caller = HybmeshCallerMatlab(exepath);
			catch
				try
					% try to load octave procedure
					ret.caller = HybmeshCallerOctave(exepath);
				catch
					error('Neither matlab nor octave caller was able to establish hybmesh server connection');
				end
			end
		end
		function free(self)
			if self.caller.connection ~= -1
				self.caller.break_connection();
				self.caller.connection = -1;
			end
		end
		function ret=apply_command(self, func, com)
			self.send_command(func, com);
			while (1)
				sig = self.wait_for_signal();
				switch (sig)
					%normal return
					case 'R'
						ret = self.read_buffer();
						return;
					%callback
					case 'B'
						self.apply_callback(self.read_buffer());
					%exception return
					case 'I'
						error('Interrupted by user');
					case 'E'
						error(self.tos_vecbyte(self.read_buffer()));
					%something went wrong
					case '0'
						error('Server stopped working');
					otherwise
						error('Invalid client instruction');
				end
			end
		end
		function ret = to_object(self, ret, str)
			ss = self.parse_vecstring(str);
			if isempty(ss)
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
		function ret = to_vecobject(self, cls, str)
			ss = self.parse_vecstring(str);
			if isempty(ss)
				ret={};
				return;
			end
			ret = {};
			ret{length(ss)} = [];
			for i=1:length(ss)
				ret{i} = self.to_object(cls, ss{i});
			end
		end
		% string -> cpp type
		function ret=to_int(~, str)
			ret = str2double(str);
		end
		function ret=to_double(~, str)
			ret = str2double(str);
		end
		function ret=to_point(~, str)
			if (strcmp(str, 'None'))
				ret = nan;
			else
				ret = str2num(str);  %#ok<ST2NM>
			end
		end
		function ret=to_grid(self, str)
			ret = self.to_object(HybmeshGrid2D, str);
		end
		function ret=to_cont(self, str)
			ret = self.to_object(HybmeshContour2D, str);
		end
		function ret=to_grid3(self, str)
			ret = self.to_object(HybmeshGrid3D, str);
		end
		function ret=to_surface(self, str)
			ret = self.to_object(HybmeshSurface3D, str);
		end
		function ret=to_veccont(self, str)
			ret = self.to_vecobject(HybmeshContour2D, str);
		end
		function ret=to_vecgrid(self, str)
			ret = self.to_vecobject(HybmeshGrid2D, str);
		end
		function ret=to_vecsurface(self, str)
			ret = self.to_vecobject(HybmeshSurface3D, str);
		end
		function ret=to_vecgrid3(self, str)
			ret = self.to_vecobject(HybmeshGrid3D, str);
		end
		function ret=to_vecint(~, str)
			ret = str2num(str);  %#ok<ST2NM>
		end
		% cpp type -> string
		function ret=tos_bool(~, val)
			if val
				ret = 'True';
			else
				ret = 'False';
			end
		end
		function ret=tos_string(self, val)
			if self.isnone(val)
				ret = 'None';
			else
				ret = strcat('''', val, '''');
			end
		end
		function ret=tos_int(~, val)
			ret =  sprintf('%i', val);
		end
		function ret=tos_double(~, val)
			ret = sprintf('%.16g', val);
		end
		function ret=tos_point(self, val)
			if (self.isnone(val))
				ret = 'None';
			else
				ret = sprintf('[%.16g, %.16g]', val);
			end
		end
		function ret=tos_point3(self, val)
			if (self.isnone(val))
				ret = 'None';
			else
				ret = sprintf('[%.16g, %.16g, %.16g]', val);
			end
		end
		function ret=tos_vecbyte(~, val)
			ret = cast(val, 'char');
		end
		function ret=tos_vecint(~, val)
			if isempty(val)
				ret = '[]';
			else
				as = sprintf('%i, ', val);
				ret = strcat('[', as, ']');
			end
		end
		function ret=tos_vecdouble(~, val)
			if isempty(val)
				ret = '[]';
			else
				as = sprintf('%.16g, ', val);
				ret = strcat('[', as, ']');
			end
		end
		function ret=tos_vecstring(~, val)
			if isempty(val)
				ret = '[]';
			else
				as = sprintf('''%s'', ', val);
				ret = strcat('[', as, ']');
			end
		end
		function ret=tos_vecpoint(self, val)
			if self.isnone(val)
				ret = 'None';
			elseif isempty(val)
				ret = '[]';
			else
				ret = sprintf('[%.16g, %.16g], ', val');
				ret = strcat('[', ret, ']');
			end
		end
		function ret=tos_object(self, val)
			if self.isnone(val)
				ret = 'None';
			else
				ret = self.tos_string(val.sid);
			end
		end
		function ret=tos_vecobject(self, val)
			if self.isnone(val)
				ret = 'None';
			elseif isempty(val)
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
		function ret=to_vec_int_double_raw(~, val)
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
		function ret=to_vecdouble_raw(~, val)
			ret = typecast(val(5:end), 'double');
		end
		function ret=to_vecint_raw(~, val)
			ret = cast(typecast(val(5:end), 'int32'), 'double');
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
			mlock;
			persistent path;
			if isempty(path)
				path = 'hybmesh.exe';  %>>$EXEPATH
			end
			if nargin
				path = newpath;
			end
			ret = path;
		end
		function ret=hybmesh_lib_path(newpath)
			mlock;
			persistent path;
			if isempty(path)
				path = 'hybmesh.so';  %>>$LIBPATH
			end
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
