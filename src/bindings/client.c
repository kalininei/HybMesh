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

#ifdef WIN32
int HybmeshClientToServer_new(const char* path, int* id){
	intptr_t h1, h2, h3, h4;
	char cmd[1000];
	HybmeshClientToServer* con;
	HANDLE sig_client2server[2];
	HANDLE sig_server2client[2];
	HANDLE data_client2server[2];
	HANDLE data_server2client[2];
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
	if (!CreatePipe(&sig_client2server[0], &sig_client2server[1], &sa, 0) || 
	    !CreatePipe(&sig_server2client[0], &sig_server2client[1], &sa, 0) ||
	    !CreatePipe(&data_client2server[0], &data_client2server[1], &sa, 0) ||
	    !CreatePipe(&data_server2client[0], &data_server2client[1], &sa, 0)){
		return 0;
	}
	/* do not inherit client handles */
	SetHandleInformation(sig_client2server[1], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(sig_server2client[0], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(data_client2server[1], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(data_server2client[0], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(sig_client2server[0], HANDLE_FLAG_INHERIT, 1);
	SetHandleInformation(sig_server2client[1], HANDLE_FLAG_INHERIT, 1);
	SetHandleInformation(data_client2server[0], HANDLE_FLAG_INHERIT, 1);
	SetHandleInformation(data_server2client[1], HANDLE_FLAG_INHERIT, 1);

	/* forking */
	ZeroMemory(&si, sizeof(si));
	ZeroMemory(&pi, sizeof(pi));
	si.cb = sizeof(si);

	h1 = (intptr_t)sig_client2server[0];
	h2 = (intptr_t)sig_server2client[1];
	h3 = (intptr_t)data_client2server[0];
	h4 = (intptr_t)data_server2client[1];

	sprintf(cmd, "%s %lld %lld %lld %lld", path, h1, h2, h3, h4);
	printf("%s\n", cmd);

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
		printf("SubProcess started\n");
	}

	/* fill con */
	con->childid = (intptr_t)pi.hProcess;
	con->sig_read = _open_osfhandle((intptr_t)sig_server2client[0], _O_APPEND|_O_RDONLY);
	con->sig_write = _open_osfhandle((intptr_t)sig_client2server[1], _O_APPEND);
	con->data_read = _open_osfhandle((intptr_t)data_server2client[0], _O_APPEND|_O_RDONLY);
	con->data_write = _open_osfhandle((intptr_t)data_client2server[1], _O_APPEND);
	con->isopen = 1;

	/* close unused handles */
	CloseHandle(pi.hThread);
	CloseHandle(sig_server2client[1]);
	CloseHandle(sig_client2server[0]);
	CloseHandle(data_server2client[1]);
	CloseHandle(data_client2server[0]);
}

#else  /* POSIX */

int HybmeshClientToServer_new(const char* path, int* id){
	HybmeshClientToServer* con;
	char s0[16], s1[16], s2[16], s3[16];
	int sig_client2server[2],  sig_server2client[2],
	    data_client2server[2], data_server2client[2];
	int childid;

	*id = get_free_index();
	if (*id == -1){
		fprintf(stderr, "no free hybmesh connections found");
		return 0;
	} else{
		con = hybmesh_connections + *id;
	}

	/* create pipes */
	if (pipe(sig_client2server)<0 || pipe(sig_server2client)<0 ||
	    pipe(data_client2server)<0 || pipe(data_server2client)<0){
		return 0;
	}

	/* forking */
	childid = fork();
	if (childid < 0) return 0;
	else if (childid == 0){
		/* child process: server */
		sprintf(s0, "%d", sig_client2server[0]);
		sprintf(s1, "%d", sig_server2client[1]);
		sprintf(s2, "%d", data_client2server[0]);
		sprintf(s3, "%d", data_server2client[1]);
		close(sig_client2server[1]);
		close(sig_server2client[0]);
		close(data_client2server[1]);
		close(data_server2client[0]);
		execl(path, path, s0, s1, s2, s3, NULL);
		return 0;
	} else {
		/* parent process: client */
		close(sig_client2server[0]);
		close(sig_server2client[1]);
		close(data_client2server[0]);
		close(data_server2client[1]);
		con->childid = childid;
		con->sig_read = sig_server2client[0];
		con->sig_write = sig_client2server[1];
		con->data_read = data_server2client[0];
		con->data_write = data_client2server[1];
		con->isopen = 1;
		return 1;
	}
}
#endif

int HybmeshClientToServer_get_signal(int id, char* sig){
	HybmeshClientToServer* con = hybmesh_connections + id;
	*sig = '0';
	PLATFORM_READ(con->sig_read, sig, 1);
	return 1;
}

int HybmeshClientToServer_get_data(int id, int* sz, char** data){
	HybmeshClientToServer* con = hybmesh_connections + id;
	PLATFORM_READ(con->data_read, (char*)(sz), 4);
	*data = (char*)malloc(*sz);
	PLATFORM_READ(con->data_read, *data, *sz);
	return 1;
}

int HybmeshClientToServer_get_data1(int id, int* sz){
	HybmeshClientToServer* con = hybmesh_connections + id;
	PLATFORM_READ(con->data_read, (char*)(sz), 4);
	return 1;
}
int HybmeshClientToServer_get_data2(int id, int sz, char* data){
	HybmeshClientToServer* con = hybmesh_connections + id;
	PLATFORM_READ(con->data_read, data, sz);
	return 1;
}

int HybmeshClientToServer_send_signal(int id, char sig){
	HybmeshClientToServer* con = hybmesh_connections + id;
	PLATFORM_WRITE(con->sig_write, &sig, 1);
	return 1;
}

int HybmeshClientToServer_send_data(int id, int sz, char* data){
	HybmeshClientToServer* con = hybmesh_connections + id;
	PLATFORM_WRITE(con->data_write, (char*)(&sz), 4);
	PLATFORM_WRITE(con->data_write, data, sz);
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
	PLATFORM_CLOSE(con->sig_read);
	PLATFORM_CLOSE(con->sig_write);
	PLATFORM_CLOSE(con->data_read);
	PLATFORM_CLOSE(con->data_write);
	con->isopen = 0;
	return 1;
}
