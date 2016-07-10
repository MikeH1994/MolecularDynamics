#ifndef _H_HELPER_H_
#define _H_HELPER_H_

#include <string>
#include <iostream>
#include <vector>
#include "GenericPlayer.h"

class Helper {
public:
	enum playerInput {
		QUIT, BLACKJACK, BLACKJACKSIM, HILO, HI, LO, YES,NO,ERR 
	};

	static playerInput getInput();
	static playerInput getInput(std::vector<playerInput> allowedInput, std::string str);
	static GenericPlayer::Moves getUserMove();
	static GenericPlayer::Moves getUserMove(std::vector <GenericPlayer::Moves> allowedInput, std::string str);
	static bool getYesNo(std::string str);
	static std::string getStr(std::string str);
	static int getInt(int min, int max, std::string str);
	static float getFloat(float min, float max, std::string str);
	static float getMin(float a, float b);
	static bool isInt(float n);
	static std::string getMoveStr(GenericPlayer::Moves move);
};

#endif