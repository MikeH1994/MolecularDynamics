#include "Helper.h"

Helper::playerInput Helper::getInput() {
	std::string input;
	std::cin >> input;

	if (input == "QUIT") {
		return QUIT;
	}
	if (input == "BLACKJACK") {
		return BLACKJACK;
	}
	if (input == "BLACKJACK_SIM") {
		return BLACKJACKSIM;
	}
	if (input == "HILO") {
		return HILO;
	}
	if (input == "HI") {
		return HI;
	}
	if (input == "LO") {
		return LO;
	}
	if (input == "YES" || input == "Y") {
		return YES;
	}
	if (input == "NO" || input == "N") {
		return NO;
	}
	return ERR;

}

Helper::playerInput Helper::getInput(std::vector<playerInput> allowedInput, std::string str) {
	bool isOk = false;
	playerInput input;

	while (isOk == false) {
		std::cout << str << std::endl;
		input = getInput();
		for (unsigned int i = 0; i < allowedInput.size(); i++) {
			// check to see if the user input is one of the allowed types
			if (input == allowedInput[i]) {
				isOk = true;
				break;
			}
		}
		if (isOk == false) {
			std::cout << "Invalid Input" << std::endl;
		}
	}
	return input;
}

GenericPlayer::Moves Helper::getUserMove() {
	//FOLD, HIT, DOUBLE, STICK, SPLIT, QUIT

	std::string input;
	std::cin >> input;

	if (input == "FOLD") {
		return GenericPlayer::FOLD;
	}
	if (input == "HIT") {
		return GenericPlayer::HIT;
	}
	if (input == "DOUBLE") {
		return GenericPlayer::DOUBLE;
	}
	if (input == "STICK") {
		return GenericPlayer::STICK;
	}
	if (input == "QUIT") {
		return GenericPlayer::QUIT;
	}
	return GenericPlayer::ERR;
}

GenericPlayer::Moves Helper::getUserMove(std::vector <GenericPlayer::Moves> allowedInput, std::string str) {
	bool isOk = false;
	GenericPlayer::Moves input = GenericPlayer::ERR;

	while (isOk == false) {
		std::cout << str << std::endl;
		input = getUserMove();
		for (unsigned int i = 0; i < allowedInput.size(); i++) {
			// check to see if the user input is one of the allowed types
			if (input == allowedInput[i]) {
				isOk = true;
				break;
			}
		}
		if (isOk == false) {
			std::cout << "Invalid Input" << std::endl;
		}
	}
	return input;
}

bool Helper::getYesNo(std::string str) {
	//return true if user enters 'yes', false for no
	std::vector<playerInput> options = { YES, NO };
	playerInput answer = getInput(options, str);
	return (answer == YES);
}

std::string Helper::getStr(std::string str) {
	std::string input;
	std::cout << str << std::endl;
	std::cin >> input;

	return input;
}

int Helper::getInt(int min, int max, std::string str) {
	int input = max + 1;
	bool validInput = false;

	std::cout << str << std::endl;
	while (!(std::cin >> input) || input < min || input>max || !isInt((float)input))
		// while input is invalid, out of range or not an integer, repeat
	{
		std::cout << "Invalid input" << std::endl;
		std::cout << str << std::endl;
		std::cin.clear();
		std::cin.ignore(INT_MAX, '\n');
	}
	return input;
}

float Helper::getFloat(float min, float max, std::string str) {
	float input = max + 1;
	bool validInput = false;

	std::cout << str << std::endl;
	while (!(std::cin >> input) || input < min || input>max)
	{
		// while input is invalid, out of range or not an integer, repeat
		std::cout << "Invalid input" << std::endl;
		std::cout << str << std::endl;
		std::cin.clear();
		std::cin.ignore(INT_MAX, '\n');
	}
	return input;
}

float Helper::getMin(float a, float b) {
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

bool Helper::isInt(float n) {
	return std::floor(std::abs(n)) == std::abs(n);
}

std::string Helper::getMoveStr(GenericPlayer::Moves move) {
	switch (move) {
	case GenericPlayer::HIT:
		return "HIT";
	case GenericPlayer::DOUBLE:
		return "DOUBLE";
	case GenericPlayer::STICK:
		return "STICK";
	case GenericPlayer::QUIT:
		return "QUIT";
	case GenericPlayer::SPLIT:
		return "SPLIT";
	case GenericPlayer::FOLD:
		return "FOLD";
	}
	return "ERR";
}