#pragma once

#include <iostream>
#include <regex>
#include <vector>

// INPUT ALPHABET
// a | !a | (a) | a&a | a|a
// (a = Inf(a), !a = Fin(a))
//
// INPUT EXAMPLE
// 0 & !1 | (1 | 2)
// Numbers represent variable indices, 0 -> variable 0.

// Debug function for printing tokens in vector<string>
void printTokens(std::vector<std::string> input){
    std::cout << "[";
    for (std::string i : input){
        std::cout << i << ", ";
    }
    std::cout << "]" << std::endl;
}

// Tokenize input string to following alphabet:
// op (!,&,|) | a | ( | )
std::vector<std::string> tokenize(std::string input){
    std::vector<std::string> result;
    std::regex inf_re("Inf");
    std::regex fin_re("Fin");
    input = std::regex_replace(input, inf_re, "");
    input = std::regex_replace(input, fin_re, "!");

    for (int i=0; i<input.size(); i++){
        std::string s = "";
        char inp = input[i];
        if (inp == ' ')
            continue;
        else if (isdigit(inp))
        {
            s += inp;
            i++;
            for (i; i<input.size(); i++){
                inp = input[i];
                if (!isdigit(inp))
                    break;
                s += inp;
            }
            i--;
        }
        else
            s += input[i];

        result.push_back(s);
    }
    return result;
}

// Helper for infix2postfix
bool isNumber(std::string s){
    for (char c : s){
        if (!isdigit(c))
            return false;
    }
    return true;
}

// Helper for infix2postfix
bool isOperator(std::string s){
    if (s == "&" || s == "|" || s == "!")
        return true;
    return false;
}

// Takes tokenized input string (in infix) and translates to postfix
std::vector<std::string> infix2postfix(std::vector<std::string> tokens){
    std::vector<std::string> opStack;
    std::vector<std::string> outputStack;

    for (std::string s : tokens){
        if (isOperator(s)){
            if (opStack.empty())
                opStack.push_back(s);
            else{
                if (opStack.back() == "("){
                    opStack.push_back(s);
                    continue;
                }
                std::string tmp = opStack.back();
                opStack.pop_back();
                outputStack.push_back(tmp);
                opStack.push_back(s);
            }
        }
        else if (isNumber(s)){
            outputStack.push_back(s);
        }
        else if (s == "("){
            opStack.push_back(s);
        }
        else{ // s == ")"
            while (true){
                if (opStack.back() == "("){
                    opStack.pop_back();
                    break;
                }
                std::string tmp = opStack.back();
                outputStack.push_back(tmp);
                opStack.pop_back();                
            }
        }
    }

    for (std::string s : opStack){
        std::string tmp = opStack.back();
        outputStack.push_back(tmp);
        opStack.pop_back();
    }

    return outputStack;
}

bool postfixEval (std::vector<std::string> postfix, std::vector<bool> colors){
    std::reverse(postfix.begin(), postfix.end());
    std::vector<bool> resStack;

    while (!postfix.empty()){
        std::string s = postfix.back();
        postfix.pop_back();

        if (isNumber(s)){
            int tmp = stoi(s);
            //std::cout << "var " << tmp << std::endl;
            resStack.push_back(colors[tmp]);
        }
        else{
            if (s == "!"){
                bool tmp = resStack.back();
                resStack.pop_back();
                //std::cout << "not " << tmp << std::endl;
                resStack.push_back(!tmp);
            }
            else if (s == "&"){
                bool tmp1 = resStack.back();
                resStack.pop_back();
                bool tmp2 = resStack.back();
                resStack.pop_back();
                //std::cout << tmp1 << " & " << tmp2 << std::endl;
                resStack.push_back(tmp1 && tmp2);
            }
            else{
                bool tmp1 = resStack.back();
                resStack.pop_back();
                bool tmp2 = resStack.back();
                resStack.pop_back();
                //std::cout << tmp1 << " & " << tmp2 << std::endl;
                resStack.push_back(tmp1 || tmp2);
            }
        }
    }
    if (resStack.size() != 1)
        std::cout << "resStack wrong size" << std::endl;
    return resStack.back();
}
