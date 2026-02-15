#include "InputTools.h"
#include <stdexcept>  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define MIN( A, B) (((A)<(B))?(A):(B))
#define MAX( A, B) (((A)>(B))?(A):(B))


string alphanumericCharacters="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.";
string lastParsedParameter,lastParsedDefinition;

extern bool lHelpUsage;

void zapspaces(char *s){
	if (s==NULL || *s=='\0') return;
	for (char *p=s; *p; p++) {
		if (!isspace(*p)) *(s++)=*p;
	}
	*s=0;
}

void ltrim(char *s){
	if (s==NULL || *s=='\0') return;
	char *p=s;
	for (; *p && isspace(*p) ; p++); // skip leading spaces
	for (; *p; p++) *(s++)=*p; // copy the remainder
	*s=0; // terminate string
}

void rtrim(char *s){
	if (s==NULL || *s=='\0') return;
	char *p=s;
	for (; *p; p++); // find end of string
	for (; p>s && isspace(*(p-1)); p--); // skip trailing spaces
	*p=0; // terminate string
}

void trim(char *s){
	if (s==NULL || *s=='\0') return;
	ltrim(s);
	rtrim(s);
}

void zapspaces(string &s){
	char *cstr = new char[s.size()+1];
	strcpy(cstr,s.c_str());
	zapspaces(cstr);
	s=string(cstr);
	delete [] cstr;
}

void ltrim(string &s){
	char *str = new char[s.size()+1];
	strcpy(str,s.c_str());
	ltrim(str);
	s=string(str);
	delete [] str;
}

void rtrim(string &s){
	char *str = new char[s.size()+1];
	strcpy(str,s.c_str());
	rtrim(str);
	s=string(str);
	delete [] str;
}

void trim(string &s){
	char *str = new char[s.size()+1];
	strcpy(str,s.c_str());
	trim(str);
	s=string(str);
	delete [] str;
}

string ltrim(const string &s){string a = s; ltrim(a); return a;}
string rtrim(const string &s){string a = s; rtrim(a); return a;}
string  trim(const string &s){string a = s;  trim(a); return a;}


bool emptyline(const char *s){
	if (s==NULL || *s=='\0') return true;
	for (const char *p=s; *p; p++) {
		if (!isspace(*p)) return false;
	}
	return true;
}

bool emptyline(const string &s){
	return emptyline(s.c_str());
}

bool isSpaceOnly(const string &s) {return emptyline(s);}
bool isWhiteSpace(const string &s) {return emptyline(s);}


bool isComment(string s,const vector<string> &comments){
	if(comments.empty()) return false;
	ltrim(s);
	for(int k=0;k<(int)comments.size();k++) {
		if(s.compare(0,comments[k].length(),comments[k])==0) return true;
	}
	return false;
}


bool isAlphaNumOnly(const string &s){
	return s.find_first_not_of(alphanumericCharacters)==std::string::npos ;
}


vector<string>	readLines(const string &fname, bool skipEmpty, bool trimLeft, bool stripComments)
{
	vector<string> lines;
	ifstream fin(fname);
	if(!fin){
		cerr<<"Error: could not open and load file: "<<fname<<endl;
		return lines;
	}
	char line[4096+1]; // 4 KB buffer
	while (fin.good()&& !fin.eof()) {
		fin.getline(line,4096);
		if(skipEmpty && emptyline(line)) continue; // skip empty lines (whitespaces only)
		if(trimLeft) ltrim(line);
		if(stripComments && isComment(line,strtokens("/ % #"))) continue; // strip out comments
		lines.push_back(string(line));
	}
	return lines;
}




vector<string> strtokens(const string &str,const char *sep){
	char * tmp = new char[str.size()+1];
	strcpy(tmp, str.c_str());
	vector<string> tokens;
	char *tok = strtok(tmp,sep);
	while(tok != NULL){
		tokens.push_back(string(tok));
		tok = strtok(NULL,sep);
	}
	delete[] tmp;
	return tokens;
}



bool isInteger(const string &str, long *val){
	long v;
	size_t idx;
	try{
		v = stol(str,&idx);
	}
	catch (const std::invalid_argument& ia) {
	// std::cerr << "Invalid argument: " << ia.what() << '\n';
		return false;
	}	
	catch (const std::out_of_range& oor) {
	// std::cerr << "Out of Range error: " << oor.what() << '\n';
		return false;
	}
	if(idx != str.size()) return false; // extra characters are present!
	if(val) *val = v;
	return true;
}

bool isNumeric(const string &str, double *val){
	double v;
	size_t idx;	
	try{
		v = stod(str,&idx);
	}

	catch (const std::invalid_argument& ia) {
	// std::cerr << "Invalid argument: " << ia.what() << '\n';
		return false;
	}	
	catch (const std::out_of_range& oor) {
	// std::cerr << "Out of Range error: " << oor.what() << '\n';
		return false;
	}	
	if(idx != str.size()) return false; // extra characters are present!
	if(val) *val = v;
	return true;
}

void educateCarriageReturn(char *str,const char *myCR){
	char CR = '\r';
	char LF = '\n';
	int nCR=0,nLF=0;
	for(char *c=str;*c;c++) {
		if (*c == CR) { if(nCR==0) { nCR++;} else {break;}}
		if (*c == LF) { if(nLF==0) { nLF++;} else {break;}}
	}
	if (nCR==0 && nLF==0) {return;} // string without newlines
	// cout<<"CR LF "<<nCR<<' '<<nLF<<endl;
	if(strlen(myCR)==1){ // 1 char => in-place replacement
		if(nCR==1 && nLF==0) { // old Mac style
			if (*myCR==CR) return; // nothing to do
			for(char *c=str;*c;c++) {if(*c==CR) *c=*myCR;} 
			return;
		}	
		if(nCR==0 && nLF==1) { // UNIX style 
			if (*myCR==LF) return; // nothing to do
			for(char *c=str;*c;c++) {if(*c==LF) *c=*myCR;} 
			return;
		}	
		if(nCR==1 && nLF==1) { // Windows style
			char *pnew = str;
			for(char *c=str;*c;c++) {
				if(*c==CR) {*pnew=*myCR; c++;} else {*pnew=*c;}
				pnew++;
			}
			*pnew='\0'; // end string
			return;
		}	
	} else {
		cerr<<"Error: this carriage return substitution not implemented yet! Request feature!"<<endl;
		cerr<<"Add new code to file "<<__FILE__<<" at line "<<__LINE__<<endl;
		exit(1);
	}
}


int getIntParam(vector<string> args,string name,int & val,int def)
{
	int ierr = getIntParamRequired(args,name,val);
	if(ierr==1) { val=def; return 0;}
	return ierr;
}

int getFloatParam(vector<string> args,string name,float & val,float def)
{
	int ierr = getFloatParamRequired(args,name,val);
	if(ierr==1) { val=def; return 0;}
	return ierr;
}

int getIntParamRequired(vector<string> args,string name,int & val)
{
	lastParsedParameter = name;
	for(int i=0;i<(int)args.size();i++){
		if(args[i]==name){
		  if(i==(int)args.size()-1) return -1; // missing value
			long lval;
			lastParsedDefinition = args[i+1];
			if(not isInteger(args[i+1],&lval)) return -2; // wrong syntax
			val = lval; // type conversion
			return 0;
		}
	}
	return 1; // not found
}

int getFloatParamRequired(vector<string> args,string name,float & val)
{
	lastParsedParameter = name;
	for(int i=0;i<(int)args.size();i++){
		if(args[i]==name){
		  if(i==(int)args.size()-1) return -1; // missing value
			lastParsedDefinition = args[i+1];
			double dval;
			if(not isNumeric(args[i+1],&dval)) return -2; // wrong syntax
			val = dval; // type conversion
			return 0;
		}
	}
	return 1; // not found
}


string pathSeparator("/");

int mkdir_c(const char* dirname) {
  string command = "mkdir -p " + string(dirname);
  return system(command.c_str());
}

int mkdir_s(string dirname){ return mkdir_c(dirname.c_str());}

int mkdir(string dirname){ return mkdir_s(dirname);}


