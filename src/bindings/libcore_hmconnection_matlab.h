#ifndef CORE_HYBMESH_CONNECTION_MATLAB_H
#define CORE_HYBMESH_CONNECTION_MATLAB_H

int require_connection(const char* path);
void send_signal(int id, char sig);
void send_data(int id, int sz, char* data);
char get_signal(int id);
int get_data1(int id);
void get_data2(int id, int sz, void* data);
void break_connection(int id);


#endif
