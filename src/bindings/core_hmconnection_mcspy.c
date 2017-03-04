#include "core_hmconnection_mcspy.h"
#include "client.h"
#include "stdlib.h"
#include "stdio.h"

static void exception(const char* str){
	fprintf(stderr, "%s", str);
	exit(EXIT_FAILURE);
}

int require_connection(const char* path){
	int ret;
	int err = HybmeshClientToServer_new(path, &ret);
	if (err == 0) exception("native require_connection() failed");
	return ret;
}

void send_signal(int id, char sig){
	int err = HybmeshClientToServer_send_signal(id, sig);
	if (err == 0) exception("native send_signal() failed");
}

void send_data(int id, int sz, char* data){
	int err = HybmeshClientToServer_send_data(id, sz, data);
	if (err == 0) exception("native send_data() failed");
}

char get_signal(int id){
	char ret;
	int err = HybmeshClientToServer_get_signal(id, &ret);
	if (err == 0) exception("native get_signal() failed");
	return ret;
}

int get_data(int id, char** data){
	int sz, err;
	err = HybmeshClientToServer_get_data(id, &sz, data);
	if (err == 0) exception("native get_data() failed");
	return sz;
}
int get_data1(int id){
	int ret, err;
	err = HybmeshClientToServer_get_data1(id, &ret);
	if (err == 0) exception("native get_data1() failed");
	return ret;
}

void get_data2(int id, int sz, char* data){
	int err = HybmeshClientToServer_get_data2(id, sz, data);
	if (err == 0) exception("native get_data2() failed");
}

void break_connection(int id){
	int err = HybmeshClientToServer_delete(id);
	if (err == 0) exception("native break_connection() failed");
}
void free_char_array(char* data){
	if (data != NULL) free(data);
}
