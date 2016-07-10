#ifndef _H_MENU_H_
#define _H_MENU_H_

#include "Helper.h"
#include "HiLo.h"
#include "Blackjack.h"

class Menu {
public:
	
	Menu();
	void run();
	void runHiLo();
	void runBlackjack();
	void runBlackjackSimulation();
	std::pair<std::vector<HumanPlayer>, std::vector<ComputerPlayer>> getPlayerVector(int nHuman, int nAI, std::vector<std::string> playerNames, float startingMoney);
	
};

#endif