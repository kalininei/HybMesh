#include "client.h"
#include "stdlib.h"
#include "stdio.h"
#include "unistd.h"

#ifdef WIN32

#include "windows.h"
#include "io.h"
#include "fcntl.h"

#define PLATFORM_READ _read
#define PLATFORM_WRITE _write
#define PLATFORM_CLOSE _close

#else

#include "sys/wait.h"
#define PLATFORM_READ read
#define PLATFORM_WRITE write
#define PLATFORM_CLOSE close

#endif

#define MAX_HYBMESH_CONNECTIONS 100
static HybmeshClientToServer hybmesh_connections[MAX_HYBMESH_CONNECTIONS];

static int get_free_index(){
	static int first_call = 1;
	int ret = -1;
	int i;
	/* allocate connections */
	if (first_call){
		for (i=0; i<MAX_HYBMESH_CONNECTIONS; ++i){
			hybmesh_connections[i].isopen = 0;
		}
		first_call = 0;
	}

	/* find first not open connector */
	for (i=0; i<MAX_HYBMESH_CONNECTIONS; ++i){
		if (!hybmesh_connections[i].isopen){
			ret = i;
			break;
		}
	}
	return ret;
}

static int find_hybmesh(const char* path, char* exepath){
	sprintf(exepath, "%s/%s", path, "hybmesh");
	if (access(exepath, X_OK) != -1) return 1;
	sprintf(exepath, "%s/%s", path, "hybmesh.py");
	if (access(exepath, X_OK) != -1) return 1;
	sprintf(exepath, "%s/%s", path, "hybmesh.exe");
	if (access(exepath, X_OK) != -1) return 1;
	return 0;
}

int HybmeshClientToServer_new(const char* path, int* id){
#ifdef WIN32
	intptr_t h1, h2, h3, h4;
	char cmd[1000], exepath[1000];
	HybmeshClientToServer* con;
	HANDLE client2server[2];
	HANDLE server2client[2];
	SECURITY_ATTRIBUTES sa;
	PROCESS_INFORMATION pi;
	STARTUPINFO si;

	*id = get_free_index();
	if (*id == -1){
		fprintf(stderr, "no free hybmesh connections found");
		return 0;
	} else{
		con = hybmesh_connections + *id;
	}

	/* create pipes */
	sa.nLength = sizeof(sa);
	sa.bInheritHandle = TRUE;
	sa.lpSecurityDescriptor = NULL;
	if (!CreatePipe(&client2server[0], &client2server[1], &sa, 0) || 
	    !CreatePipe(&server2client[0], &server2client[1], &sa, 0)){
		return 0;
	}
	/* do not inherit client handles */
	SetHandleInformation(client2server[1], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(server2client[0], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(client2server[0], HANDLE_FLAG_INHERIT, 1);
	SetHandleInformation(server2client[1], HANDLE_FLAG_INHERIT, 1);

	/* forking */
	ZeroMemory(&si, sizeof(si));
	ZeroMemory(&pi, sizeof(pi));
	si.cb = sizeof(si);

	h1 = (intptr_t)client2server[0];
	h2 = (intptr_t)server2client[1];

	if (find_hybmesh(path, exepath) == 0){
		fprintf(stderr, "hybmesh executable is not found at %s\n", path);
		return 0;
	}
	sprintf(cmd, "%s -px %lld %lld", exepath, h1, h2);

	if(!CreateProcess( NULL,/* No module name (use command line) */
		cmd,            /* Command line*/
		NULL,           /* Process handle not inheritable */
		NULL,           /* Thread handle not inheritable */
		TRUE,           /* Set handle inheritance to TRUE */
		0,              /* No creation flags */
		NULL,           /* Use parent's environment block */
		NULL,           /* Use parent's starting directory */
		&si,            /* Pointer to STARTUPINFO structure */
		&pi )           /* Pointer to PROCESS_INFORMATION structure */
	){
		fprintf(stderr, "CreateProcess failed (%d).\n", GetLastError() );
		return 0;
	} else{
	}

	/* fill con */
	con->childid = (intptr_t)pi.hProcess;
	con->pipe_read = _open_osfhandle((intptr_t)server2client[0], _O_APPEND|_O_RDONLY);
	con->pipe_write = _open_osfhandle((intptr_t)client2server[1], _O_APPEND);
	con->isopen = 1;

	/* close unused handles */
	CloseHandle(pi.hThread);
	CloseHandle(server2client[1]);
	CloseHandle(client2server[0]);

#else  /* POSIX */

	HybmeshClientToServer* con;
	char s0[16], s1[16], exepath[1000];
	int client2server[2],  server2client[2];
	int childid;

	*id = get_free_index();
	if (*id == -1){
		fprintf(stderr, "no free hybmesh connections found");
		return 0;
	} else{
		con = hybmesh_connections + *id;
	}

	/* create pipes */
	if (pipe(client2server)<0 || pipe(server2client)<0)
		return 0;

	/* forking */
	childid = fork();
	if (childid < 0) return 0;
	else if (childid == 0){
		/* child process: server */
		sprintf(s0, "%d", client2server[0]);
		sprintf(s1, "%d", server2client[1]);
		close(client2server[1]);
		close(server2client[0]);
		if (find_hybmesh(path, exepath) == 0){
			fprintf(stderr, "hybmesh executable is not found at %s\n", path);
			return 0;
		}
		execl(exepath, exepath, "-px", s0, s1, NULL);
		return 0;
	} else {
		/* parent process: client */
		close(client2server[0]);
		close(server2client[1]);
		con->childid = childid;
		con->pipe_read = server2client[0];
		con->pipe_write = client2server[1];
		con->isopen = 1;
		return 1;
	}
#endif
}

static int read_nbytes(int fd, void* buf, size_t sz){
	int rd;

	while (1){
		rd = PLATFORM_READ(fd, buf, sz);
		if (rd < 0 || rd > sz){
			fprintf(stderr, "native read error\n");
			return 0;
		} else if (sz == rd) return 1;
		sz -= rd;
		buf += rd;
	}
}

int HybmeshClientToServer_get_signal(int id, char* sig){
	HybmeshClientToServer* con = hybmesh_connections + id;
	*sig = '0';
	if (read_nbytes(con->pipe_read, sig, 1) == 0) return 0;
	return 1;
}

int HybmeshClientToServer_get_data(int id, int* sz, char** data){
	HybmeshClientToServer* con = hybmesh_connections + id;
	if (read_nbytes(con->pipe_read, sz, 4) == 0) return 0;
	*data = (char*)malloc(*sz);
	if (read_nbytes(con->pipe_read, *data, *sz) == 0) return 0;
	return 1;
}

int HybmeshClientToServer_get_data1(int id, int* sz){
	HybmeshClientToServer* con = hybmesh_connections + id;
	if (read_nbytes(con->pipe_read, sz, 4) == 0) return 0;
	return 1;
}
int HybmeshClientToServer_get_data2(int id, int sz, char* data){
	HybmeshClientToServer* con = hybmesh_connections + id;
	if (read_nbytes(con->pipe_read, data, sz) == 0) return 0;
	return 1;
}

int HybmeshClientToServer_send_signal(int id, char sig){
	HybmeshClientToServer* con = hybmesh_connections + id;
	PLATFORM_WRITE(con->pipe_write, &sig, 1);
	return 1;
}

int HybmeshClientToServer_send_data(int id, int sz, char* data){
	HybmeshClientToServer* con = hybmesh_connections + id;
	PLATFORM_WRITE(con->pipe_write, &sz, 4);
	PLATFORM_WRITE(con->pipe_write, data, sz);
	return 1;
}

int HybmeshClientToServer_delete(int id){
	HybmeshClientToServer* con = hybmesh_connections + id;
	HybmeshClientToServer_send_signal(id, 'Q');
#ifdef WIN32
	WaitForSingleObject((HANDLE)con->childid, INFINITE);
#else
	waitpid(con->childid, NULL, 0);
#endif
	PLATFORM_CLOSE(con->pipe_read);
	PLATFORM_CLOSE(con->pipe_write);
	con->isopen = 0;
	return 1;
}
