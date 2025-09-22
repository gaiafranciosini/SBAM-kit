#pragma once

#include <cstring>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
using namespace std;

extern string alphanumericCharacters;
extern string lastParsedParameter,lastParsedDefinition;

void	zapspaces(string &s);

vector<string> readLines(const string &fname, bool skipEmpty=true, bool trimLeft=true, bool stripComments=true);

vector<string> strtokens(const string &str,const char *sep = " \t");

string ltrim(const string &s);
string rtrim(const string &s);
string trim(const string &s);

void ltrim(string &s);
void rtrim(string &s);
void  trim(string &s);

bool emptyline(const string &s);

bool isInteger(const string &s, long *val=NULL);
bool isNumeric(const string &s, double *val=NULL);

bool isSpaceOnly(const string &s);
bool isWhiteSpace(const string &s);
bool isComment(string s,const vector<string> &comments);
bool isAlphaNumOnly(const string &s);

void educateCarriageReturn(char *str,const char *myCR);


int getIntParam(vector<string> args,string name,int & val,int def);
int getFloatParam(vector<string> args,string name,float & val,float def);

int getIntParamRequired(vector<string> args,string name,int & val);
int getFloatParamRequired(vector<string> args,string name,float & val);

int mkdir(string dirname);

