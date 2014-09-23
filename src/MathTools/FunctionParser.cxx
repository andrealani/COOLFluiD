// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

//===============================
// Function parser v2.22 by Warp
//===============================

//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <new>
#include <cmath>

#include "MathTools/FunctionParser.hh"

#ifndef CF_HAVE_CUDA
#ifdef CF_HAVE_BOOST_ERFC
  #include <boost/math/special_functions/erf.hpp>
#endif
#endif

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

    const char* const FuncNames[]=
    {
        "abs","exp","ceil","floor","log","sqrt","int",
        "sinh","cosh","tanh","sin","cos","tan",
        "arctan",
#if   ( defined(CF_HAVE_BOOST_ERFC) || defined(CF_HAVE_MATH_ERFC) )
        "erfc",
#endif
#ifdef CF_HAVE_MATH_ASINH
        "asinh",
#endif
#ifdef CF_HAVE_MATH_ACOSH
        "acosh",
#endif
#ifdef CF_HAVE_MATH_ATANH
        "atanh",
#endif
        "asin","acos","atan",
        "min", "max", "if",
#ifdef CF_ENABLE_FUNCTIONPARSER_EVAL
        "eval",
#endif
        0
    };

    const unsigned FuncParams[]=
    {
        1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1,
        2,
#if   ( defined (CF_HAVE_BOOST_ERFC) || defined(CF_HAVE_MATH_ERFC) )
        1,
#endif
#ifdef CF_HAVE_MATH_ASINH
        1,
#endif
#ifdef CF_HAVE_MATH_ACOSH
        1,
#endif
#ifdef CF_HAVE_MATH_ATANH
        1,
#endif
        1, 1, 1,
        2, 2, 0,
#ifdef CF_ENABLE_FUNCTIONPARSER_EVAL
        0
#endif
    };

    enum OPCODE
    {
        cImmed, cJump,
        cNeg, cAdd, cSub, cMul, cDiv, cMod, cPow,
        cEqual, cLess, cGreater, cAnd, cOr,
        cAbs, cExp, cCeil, cFloor, cLog, cSqrt, cInt,
        cSinh, cCosh, cTanh, cSin, cCos, cTan,
        cArctan,
#if   ( defined(CF_HAVE_BOOST_ERFC) || defined(CF_HAVE_MATH_ERFC) )
        cErfc,
#endif
#ifdef CF_HAVE_MATH_ASINH
        cAsinh,
#endif
#ifdef CF_HAVE_MATH_ACOSH
        cAcosh,
#endif
#ifdef CF_HAVE_MATH_ATANH
        cAtanh,
#endif
        cAsin, cAcos, cAtan,
        cMin, cMax, cIf,
#ifdef CF_ENABLE_FUNCTIONPARSER_EVAL
        cEval,
#endif
        VarBegin
    };

    struct FuncDefinition
    {
        unsigned opcode;
        unsigned params;

        inline FuncDefinition(unsigned o, unsigned p): opcode(o), params(p) {}
        inline FuncDefinition(): opcode(0), params(1) {}
    };

    typedef map<string, FuncDefinition> FuncMap_t;
    FuncMap_t Functions;

    template<typename MapType>
    inline typename MapType::const_iterator
    FindInMap(const MapType& Map, const char* F)
    {
        unsigned ind = 0;
        while(isalpha(F[ind])) ++ind;
        if(ind)
        {
            string name(F, ind);
            return Map.find(name);
        }
        return Map.end();
    }

    inline FuncMap_t::const_iterator FindFunction(const char* F)
    {
        return FindInMap(Functions, F);
    }

//////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------
// Constructors and destructors
//---------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
FunctionParser::FunctionParser() : ParseErrorType(-1), EvalErrorType(0)
{
    // Initialize function name map
    if(Functions.size() == 0)
    {
        for(unsigned fInd = 0; FuncNames[fInd]; ++fInd)
        {
            Functions[FuncNames[fInd]] =
                FuncDefinition(cAbs+fInd, FuncParams[fInd]);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

// Copy constructor (only for private use)
FunctionParser::FunctionParser(const FunctionParser& cpy):
    varAmount(cpy.varAmount),
    EvalErrorType(cpy.EvalErrorType),
    Comp(cpy.Comp)
{
}

//////////////////////////////////////////////////////////////////////////////

FunctionParser::~FunctionParser()
{
}

//////////////////////////////////////////////////////////////////////////////

FunctionParser::CompiledCode::CompiledCode():
    ByteCode(0), ByteCodeSize(0),
    Immed(0), ImmedSize(0),
    Stack(0), StackSize(0),
    thisIsACopy(false)
{
}

//////////////////////////////////////////////////////////////////////////////

FunctionParser::CompiledCode::CompiledCode(const CompiledCode& cpy):
    ByteCode(cpy.ByteCode), ByteCodeSize(cpy.ByteCodeSize),
    Immed(cpy.Immed), ImmedSize(cpy.ImmedSize),
    Stack(new double[cpy.StackSize]), StackSize(cpy.StackSize),
    thisIsACopy(true)
{
}

//////////////////////////////////////////////////////////////////////////////

FunctionParser::CompiledCode::~CompiledCode()
{
    if(!thisIsACopy && ByteCode) { delete[] ByteCode; ByteCode=0; }
    if(!thisIsACopy && Immed) { delete[] Immed; Immed=0; }
    if(Stack) { delete[] Stack; Stack=0; }
}

//////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------
// Function parsing
//---------------------------------------------------------------------------

    const char* ParseErrorMessage[]=
    {
        "Syntax error",                             // 0
        "Mismatched parenthesis",                   // 1
        "Missing ')'",                              // 2
        "Empty parentheses",                        // 3
        "Syntax error. Operator expected",          // 4
        "Not enough memory",                        // 5
        "An unexpected error ocurred. Please make a full bug report "
        "to warp@iki.fi",                           // 6
        "Syntax error in parameter 'Vars' given to "
        "FunctionParser::Parse()",                  // 7
        "Illegal number of parameters to function"  // 8
    };

    // Return index to original string when the index is to the function
    // with no whitespaces
    inline int RealIndexPos(const string& s,int Pos)
    {
        int i, Ind=0;
        for(i=0; i<Pos; i++,Ind++)
            while(s[Ind]==' ') Ind++;
        while(s[Ind]==' ') Ind++;
        return Ind;
    }

    typedef map<string, unsigned> NameMap_t;

    // Parse variables
    bool ParseVars(const string& Vars, NameMap_t& dest)
    {
        unsigned varNumber = VarBegin;
        unsigned ind1 = 0, ind2;

        while(ind1 < Vars.size())
        {
            for(ind2=ind1; ind2<Vars.size() && Vars[ind2]!=','; ++ind2) // loops until next comma, points ind2 to end of var
                if(!isalpha(Vars[ind2])) return false;
            if(ind2 == ind1) return false; // detects consecutive commas 
            const string varName = Vars.substr(ind1, ind2-ind1);

          
            /// detect duplicated variable names
            NameMap_t::const_iterator iter = dest.find(varName);
            if(iter != dest.end()) return false;
            dest[varName] = varNumber++;

            ind1 = ind2+1;
        }
        return true;
    }

//////////////////////////////////////////////////////////////////////////////

// Main parsing function
int FunctionParser::Parse(const std::string& Function, const std::string& Vars)
{
    Variables.clear();

    if(!ParseVars(Vars, Variables))
    {
        ParseErrorType = 7;
        return Function.size();
    }
    varAmount = Variables.size(); // this is for Eval()

    string tmp;
    tmp.reserve(Function.size());
    for(unsigned i=0; i<Function.size(); ++i)
        if(!isspace(Function[i])) tmp += Function[i];
    const char* Func = tmp.c_str();

    ParseErrorType = -1;

    int Result = CheckSyntax(Func);
    if(Result>=0)
    {
        return RealIndexPos(Function, Result);
    }

    if(!Compile(Func)) return Function.size();

    Variables.clear();

    ParseErrorType = -1;
    return -1;
}

//////////////////////////////////////////////////////////////////////////////

    // Is given char an operator?
    inline bool IsOperator(int c)
    {
        return strchr("+-*/%^=<>&|,",c)!=NULL;
    }

//////////////////////////////////////////////////////////////////////////////

inline FunctionParser::VarMap_t::const_iterator
FunctionParser::FindVariable(const char* F)
{
    return FindInMap(Variables, F);
}

//////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------
// Check function string syntax
// ----------------------------
int FunctionParser::CheckSyntax(const char* Function)
{
    int Ind=0,ParenthCnt=0,c;
    char* Ptr;

    while(true)
    {
        c=Function[Ind];

// Check for valid operand (must appear)

        // Check for leading -
        if(c=='-') c=Function[++Ind];
        if(c==0) { ParseErrorType=0; return Ind; }

        // Check for math function
        FuncMap_t::const_iterator fIter = FindFunction(&Function[Ind]);
        if(fIter != Functions.end())
        {
            Ind += fIter->first.size();
            c = Function[Ind];
            if(c!='(') { ParseErrorType=0; return Ind; }
        }

        // Check for opening parenthesis
        if(c=='(') { ParenthCnt++; Ind++; continue; }

        // Check for number
        if(isdigit(c) || (c=='.' && isdigit(Function[Ind+1])))
        {
            strtod(&Function[Ind], &Ptr);
            Ind += int(Ptr-&Function[Ind]);
            c = Function[Ind];
        }
        else
        { // Check for variable
            VarMap_t::const_iterator vIter = FindVariable(&Function[Ind]);
            if(vIter == Variables.end()) { ParseErrorType=0; return Ind; }
            Ind += vIter->first.size();
            c = Function[Ind];
        }

        // Check for closing parenthesis
        while(c==')')
        {
            ParenthCnt--;
            if(ParenthCnt<0) { ParseErrorType=1; return Ind; }
            if(Function[Ind-1]=='(') { ParseErrorType=3; return Ind; }
            c=Function[++Ind];
        }

// If we get here, we have a legal operand and now a legal operator or
// end of string must follow

        // Check for EOS
        if(c==0) break; // The only way to end the checking loop without error
        // Check for operator
        if(!IsOperator(c)) { ParseErrorType=4; return Ind; }

// If we get here, we have an operand and an operator; the next loop will
// check for another operand (must appear)
        ++Ind;
    } // while

    // Check that all opened parentheses are also closed
    if(ParenthCnt>0) { ParseErrorType=2; return Ind; }

// The string is ok
    ParseErrorType=-1;
    return -1;
}

//////////////////////////////////////////////////////////////////////////////

// Compile function string to bytecode
// -----------------------------------
bool FunctionParser::Compile(const char* Function)
{
    if(Comp.ByteCode) { delete[] Comp.ByteCode; Comp.ByteCode=0; }
    if(Comp.Immed) { delete[] Comp.Immed; Comp.Immed=0; }
    if(Comp.Stack) { delete[] Comp.Stack; Comp.Stack=0; }

// Compile to nowhere to get the size of the bytecode
    Comp.ByteCodeSize=Comp.ImmedSize=Comp.StackSize=Comp.StackPtr=0;
    CompileExpression(Function, 0);
    if(ParseErrorType >= 0) return false;

    try
    {
        if(Comp.ByteCodeSize)
            Comp.ByteCode = new unsigned[Comp.ByteCodeSize];
        if(Comp.ImmedSize)
            Comp.Immed = new double[Comp.ImmedSize];
        if(Comp.StackSize)
            Comp.Stack = new double[Comp.StackSize];
    }
    catch(std::bad_alloc)
    {
        ParseErrorType=5; return false;
        // Already allocated memory blocks are freed in the destructor or when
        // Parse() is called again, so no need to free them here
    }

// Compile
    Comp.ByteCodeSize=Comp.ImmedSize=Comp.StackSize=Comp.StackPtr=0;
    CompileExpression(Function, 0);

    return ParseErrorType < 0;
}

//////////////////////////////////////////////////////////////////////////////

inline void FunctionParser::AddCompiledByte(unsigned c)
{
    if(Comp.ByteCode)
        Comp.ByteCode[Comp.ByteCodeSize] = c;
    ++Comp.ByteCodeSize;
}

//////////////////////////////////////////////////////////////////////////////

// Compile if()
int FunctionParser::CompileIf(const char* F, int ind)
{
    int ind2 = CompileExpression(F, ind, true); // condition
    if(F[ind2] != ',') { ParseErrorType=8; return ind2; }
    AddCompiledByte(cIf);
    unsigned curByteCodeSize = Comp.ByteCodeSize;
    AddCompiledByte(0); // Jump index; to be set later
    AddCompiledByte(0); // Immed jump index; to be set later

    --Comp.StackPtr;

    ind2 = CompileExpression(F, ind2+1, true); // then
    if(F[ind2] != ',') { ParseErrorType=8; return ind2; }
    AddCompiledByte(cJump);
    unsigned curByteCodeSize2 = Comp.ByteCodeSize;
    unsigned curImmedSize2 = Comp.ImmedSize;
    AddCompiledByte(0); // Jump index; to be set later
    AddCompiledByte(0); // Immed jump index; to be set later

    --Comp.StackPtr;

    ind2 = CompileExpression(F, ind2+1, true); // else
    if(F[ind2] != ')') { ParseErrorType=8; return ind2; }

    if(Comp.ByteCode) // Set jump indices
    {
        Comp.ByteCode[curByteCodeSize] = curByteCodeSize2+1;
        Comp.ByteCode[curByteCodeSize+1] = curImmedSize2;
        Comp.ByteCode[curByteCodeSize2] = Comp.ByteCodeSize-1;
        Comp.ByteCode[curByteCodeSize2+1] = Comp.ImmedSize;
    }

    return ind2+1;
}

//////////////////////////////////////////////////////////////////////////////

// Compiles element
int FunctionParser::CompileElement(const char* F, int ind)
{
    char c = F[ind];

    if(c == '(')
    {
        return CompileExpression(F, ind+1) + 1;
    }
    else if(c == '-')
    {
        int ind2 = CompileElement(F, ind+1);
        AddCompiledByte(cNeg);
        return ind2;
    }

    if(isdigit(c) || c=='.') // Number
    {
        const char* startPtr = &F[ind];
        char* endPtr;
        double val = strtod(startPtr, &endPtr);
        if(Comp.Immed)
            Comp.Immed[Comp.ImmedSize] = val;
        ++Comp.ImmedSize;
        AddCompiledByte(cImmed);
        ++Comp.StackPtr; if(Comp.StackPtr>Comp.StackSize) Comp.StackSize++;
        return ind+(endPtr-startPtr);
    }

    if(isalpha(c)) // Function or variable
    {
        int ind2 = ind+1;
        while(isalpha(F[ind2])) ++ind2;
        string name(F+ind, ind2-ind);

        FuncMap_t::const_iterator fIter = Functions.find(name);
        if(fIter != Functions.end()) // is function
        {
            if(fIter->first == "if") // "if" is a special case
            {
                return CompileIf(F, ind2+1);
            }

            unsigned curStackPtr = Comp.StackPtr;
            ind2 = CompileExpression(F, ind2+1);

#ifdef CF_ENABLE_FUNCTIONPARSER_EVAL
            unsigned requiredParams =
                fIter->first=="eval" ? Variables.size() : fIter->second.params;
#else
            unsigned requiredParams = fIter->second.params;
#endif
            if(Comp.StackPtr != curStackPtr+requiredParams)
            { ParseErrorType=8; return ind; }

            AddCompiledByte(fIter->second.opcode);
            Comp.StackPtr -= fIter->second.params - 1;
            return ind2+1;
        }

        VarMap_t::const_iterator vIter = Variables.find(name);
        if(vIter != Variables.end()) // is variable
        {
            AddCompiledByte(vIter->second);
            ++Comp.StackPtr; if(Comp.StackPtr>Comp.StackSize) Comp.StackSize++;
            return ind + vIter->first.size();
        }
    }

    ParseErrorType = 6;
    return ind;
}

//////////////////////////////////////////////////////////////////////////////

// Compiles '^'
int FunctionParser::CompilePow(const char* F, int ind)
{
    int ind2 = CompileElement(F, ind);

    while(F[ind2] == '^')
    {
        ind2 = CompileElement(F, ind2+1);
        AddCompiledByte(cPow);
        --Comp.StackPtr;
    }

    return ind2;
}

//////////////////////////////////////////////////////////////////////////////

// Compiles '*', '/' and '%'
int FunctionParser::CompileMult(const char* F, int ind)
{
    int ind2 = CompilePow(F, ind);
    char op;

    while((op = F[ind2]) == '*' || op == '/' || op == '%')
    {
        ind2 = CompilePow(F, ind2+1);
        switch(op)
        {
          case '*': AddCompiledByte(cMul); break;
          case '/': AddCompiledByte(cDiv); break;
          case '%': AddCompiledByte(cMod); break;
        }
        --Comp.StackPtr;
    }

    return ind2;
}

//////////////////////////////////////////////////////////////////////////////

// Compiles '+' and '-'
int FunctionParser::CompileAddition(const char* F, int ind)
{
    int ind2 = CompileMult(F, ind);
    char op;

    while((op = F[ind2]) == '+' || op == '-')
    {
        ind2 = CompileMult(F, ind2+1);
        AddCompiledByte(op=='+' ? cAdd : cSub);
        --Comp.StackPtr;
    }

    return ind2;
}

//////////////////////////////////////////////////////////////////////////////

// Compiles '=', '<' and '>'
int FunctionParser::CompileComparison(const char* F, int ind)
{
    int ind2 = CompileAddition(F, ind);
    char op;

    while((op = F[ind2]) == '=' || op == '<' || op == '>')
    {
        ind2 = CompileAddition(F, ind2+1);
        switch(op)
        {
          case '=': AddCompiledByte(cEqual); break;
          case '<': AddCompiledByte(cLess); break;
          case '>': AddCompiledByte(cGreater); break;
        }
        --Comp.StackPtr;
    }

    return ind2;
}

//////////////////////////////////////////////////////////////////////////////

// Compiles '&'
int FunctionParser::CompileAnd(const char* F, int ind)
{
    int ind2 = CompileComparison(F, ind);

    while(F[ind2] == '&')
    {
        ind2 = CompileComparison(F, ind2+1);
        AddCompiledByte(cAnd);
        --Comp.StackPtr;
    }

    return ind2;
}

//////////////////////////////////////////////////////////////////////////////

// Compiles '|'
int FunctionParser::CompileOr(const char* F, int ind)
{
    int ind2 = CompileAnd(F, ind);

    while(F[ind2] == '|')
    {
        ind2 = CompileAnd(F, ind2+1);
        AddCompiledByte(cOr);
        --Comp.StackPtr;
    }

    return ind2;
}

//////////////////////////////////////////////////////////////////////////////

// Compiles ','
int FunctionParser::CompileExpression(const char* F, int ind, bool stopAtComma)
{
    int ind2 = CompileOr(F, ind);

    if(stopAtComma) return ind2;

    while(F[ind2] == ',')
    {
        ind2 = CompileOr(F, ind2+1);
    }

    return ind2;
}

//////////////////////////////////////////////////////////////////////////////


// Return parse error message
// --------------------------
const char* FunctionParser::ErrorMsg(void) const
{
    if(ParseErrorType>=0) return ParseErrorMessage[ParseErrorType];
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------
// Function evaluation
//---------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////

    inline int doubleToInt(double d)
    {
        return d<0 ? -int((-d)+.5) : int(d+.5);
    }

    inline double Min(double d1, double d2)
    {
        return d1<d2 ? d1 : d2;
    }
    inline double Max(double d1, double d2)
    {
        return d1>d2 ? d1 : d2;
    }

//////////////////////////////////////////////////////////////////////////////

//double FunctionParser::Eval(const double* Vars)
double FunctionParser::Eval(const RealVector& Vars)
{
    unsigned IP, DP=0;
    int SP=-1;

    for(IP=0; IP<Comp.ByteCodeSize; IP++)
    {
        switch(Comp.ByteCode[IP])
        {
          case cImmed: Comp.Stack[++SP]=Comp.Immed[DP++]; break;

          case  cJump: DP = Comp.ByteCode[IP+2];
                       IP = Comp.ByteCode[IP+1];
                       break;

          case   cNeg: Comp.Stack[SP]=-Comp.Stack[SP]; break;
          case   cAdd: Comp.Stack[SP-1]+=Comp.Stack[SP]; SP--; break;
          case   cSub: Comp.Stack[SP-1]-=Comp.Stack[SP]; SP--; break;
          case   cMul: Comp.Stack[SP-1]*=Comp.Stack[SP]; SP--; break;
          case   cDiv: if(Comp.Stack[SP]==0) { EvalErrorType=1; return 0; }
                       Comp.Stack[SP-1]/=Comp.Stack[SP]; SP--; break;
          case   cMod: if(Comp.Stack[SP]==0) { EvalErrorType=1; return 0; }
                       Comp.Stack[SP-1]=fmod(Comp.Stack[SP-1],Comp.Stack[SP]);
                       SP--; break;
          case   cPow: Comp.Stack[SP-1]=pow(Comp.Stack[SP-1],Comp.Stack[SP]);
                       SP--; break;

          case cEqual: Comp.Stack[SP-1] = (Comp.Stack[SP-1]==Comp.Stack[SP]);
                       SP--; break;
          case  cLess: Comp.Stack[SP-1] = (Comp.Stack[SP-1]<Comp.Stack[SP]);
                       SP--; break;
          case cGreater: Comp.Stack[SP-1] = (Comp.Stack[SP-1]>Comp.Stack[SP]);
                         SP--; break;
          case   cAnd: Comp.Stack[SP-1] =
                           (doubleToInt(Comp.Stack[SP-1]) &&
                            doubleToInt(Comp.Stack[SP]));
                       SP--; break;
          case    cOr: Comp.Stack[SP-1] =
                           (doubleToInt(Comp.Stack[SP-1]) ||
                            doubleToInt(Comp.Stack[SP]));
                       SP--; break;

          case   cAbs: Comp.Stack[SP]=std::abs(Comp.Stack[SP]); break;
	case   cExp: Comp.Stack[SP]=std::exp(Comp.Stack[SP]); break;
          case  cCeil: Comp.Stack[SP]=ceil(Comp.Stack[SP]); break;
          case cFloor: Comp.Stack[SP]=floor(Comp.Stack[SP]); break;
          case   cLog: if(Comp.Stack[SP]<=0) { EvalErrorType=3; return 0; }
	    Comp.Stack[SP]=std::log(Comp.Stack[SP]); break;
          case  cSqrt: if(Comp.Stack[SP]<0) { EvalErrorType=2; return 0; }
	    Comp.Stack[SP]=std::sqrt(Comp.Stack[SP]); break;
          case   cInt: Comp.Stack[SP]=doubleToInt(Comp.Stack[SP]); break;
	case  cSinh: Comp.Stack[SP]=std::sinh(Comp.Stack[SP]); break;
	case  cCosh: Comp.Stack[SP]=std::cosh(Comp.Stack[SP]); break;
	case  cTanh: Comp.Stack[SP]=std::tanh(Comp.Stack[SP]); break;
	case   cSin: Comp.Stack[SP]=std::sin(Comp.Stack[SP]); break;
	case   cCos: Comp.Stack[SP]=std::cos(Comp.Stack[SP]); break;
	case   cTan: Comp.Stack[SP]=std::tan(Comp.Stack[SP]); break;
          case cArctan: Comp.Stack[SP-1]=atan2(Comp.Stack[SP-1],Comp.Stack[SP]);
                       SP--; break;

#ifndef CF_HAVE_CUDA 
  #ifdef CF_HAVE_BOOST_ERFC
          case  cErfc: Comp.Stack[SP]=boost::math::erfc(Comp.Stack[SP]); break;
  #else
     #ifdef CF_HAVE_MATH_ERFC
          case  cErfc: Comp.Stack[SP]=erfc(Comp.Stack[SP]); break;
     #endif
  #endif
#else
     #ifdef CF_HAVE_MATH_ERFC
          case  cErfc: Comp.Stack[SP]=erfc(Comp.Stack[SP]); break;
     #endif
#endif

#ifdef CF_HAVE_MATH_ASINH
          case cAsinh: Comp.Stack[SP]=asinh(Comp.Stack[SP]); break;
#endif
#ifdef CF_HAVE_MATH_ACOSH
          case cAcosh: Comp.Stack[SP]=acosh(Comp.Stack[SP]); break;
#endif
#ifdef CF_HAVE_MATH_ATANH
          case cAtanh: Comp.Stack[SP]=atanh(Comp.Stack[SP]); break;
#endif
          case  cAsin: if(Comp.Stack[SP]<-1 || Comp.Stack[SP]>1)
                       { EvalErrorType=4; return 0; }
	    Comp.Stack[SP]=std::asin(Comp.Stack[SP]); break;
          case  cAcos: if(Comp.Stack[SP]<-1 || Comp.Stack[SP]>1)
                       { EvalErrorType=4; return 0; }
	    Comp.Stack[SP]=std::acos(Comp.Stack[SP]); break;
	case  cAtan: Comp.Stack[SP]=std::atan(Comp.Stack[SP]); break;

          case   cMin: Comp.Stack[SP-1]=Min(Comp.Stack[SP-1],Comp.Stack[SP]);
                       SP--; break;
          case   cMax: Comp.Stack[SP-1]=Max(Comp.Stack[SP-1],Comp.Stack[SP]);
                       SP--; break;
          case    cIf:
              {
                  unsigned jumpAddr = Comp.ByteCode[++IP];
                  unsigned immedAddr = Comp.ByteCode[++IP];
                  if(doubleToInt(Comp.Stack[SP]) == 0)
                  {
                      IP = jumpAddr;
                      DP = immedAddr;
                  }
                  SP--; break;
              }

#ifdef CF_ENABLE_FUNCTIONPARSER_EVAL
          case  cEval:
              {
                  FunctionParser fpcopy(*this);
                  double retVal = fpcopy.Eval(&Comp.Stack[SP-varAmount+1]);
                  SP -= varAmount-1;
                  Comp.Stack[SP] = retVal;
                  break;
              }
#endif

          default:
              Comp.Stack[++SP]=Vars[Comp.ByteCode[IP]-VarBegin];
        }
    }

    EvalErrorType=0;
    return Comp.Stack[SP];
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

