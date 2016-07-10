#include "Menu.h"

Menu::Menu() {
	while (true) {
		run();
	}
}

void Menu::run() {
	std::cout << std::endl;
	Helper::playerInput input = Helper::ERR;
	std::vector<Helper::playerInput> allowedOptions = { Helper::BLACKJACK, Helper::BLACKJACKSIM, Helper::HILO };
	while (input == Helper::ERR) {
		input = Helper::getInput(allowedOptions, "Select a mode to run- 'BLACKJACK', 'BLACKJACK_SIM' or 'HILO'");
	}
	switch (input) {
	case Helper::BLACKJACK:
		runBlackjack();
		break;
	case Helper::BLACKJACKSIM:
		runBlackjackSimulation();
		break;
	case Helper::HILO:
		runHiLo();
		break;
	}
}

void Menu::runHiLo() {
	HiLo* hilo = new HiLo(2);
}

void Menu::runBlackjack() {
	std::cout << "======================Settings======================" << std::endl;
	float startingMoney = 20;
	int nHuman = Helper::getInt(0, 2, "Enter the number of human players (max 2)");
	std::vector<std::string> nameVec(nHuman);
	for (int i = 0; i < nHuman; i++) {
		nameVec[i] = Helper::getStr("Enter name for player " + std::to_string(i + 1));
	}
	int nMin = 0;
	if (nHuman == 0) {
		nMin = 1;
	}
	int nAI = Helper::getInt(0, 5, "Enter the number of AI players (max 5)");
	int nDecks = Helper::getInt(2, 8, "Enter the number of decks (min 2, max 8)");

	std::cout << "\nHuman Players: " << nHuman << "\nAI Players: " << nAI << "\nNumber of decks: " << nDecks << std::endl;

	bool _continue = Helper::getYesNo("Continue with options? (YES/NO)");
	if (_continue) {
		// if the user is happy with options, continue
		std::pair<std::vector<HumanPlayer>, std::vector<ComputerPlayer>> pair = getPlayerVector(nHuman, nAI, nameVec, startingMoney);
		Blackjack* blackjack = new Blackjack(pair.first, pair.second, nDecks, true, true);
		blackjack->run();
	}
	else {
		// if the user chooses 'no'
		runBlackjack();
	}

}

void Menu::runBlackjackSimulation() {
	std::cout << "======================Settings======================" << std::endl;
	int nTimes = Helper::getInt(1, (int) 1E7, "Enter the number of hands to run");
	int nDecks = Helper::getInt(1, 8, "Enter the number of decks (min 1, max 8)");
	float winPay = Helper::getFloat(0, 1000, "Enter return ratio for a win (not including initial bet)");
	float blackjackPay = Helper::getFloat(0, 1000, "Enter return ratio for a blackjack (not including initial bet)");

	std::cout << "Number of hands: " << nTimes << "\nNumber of decks: " << nDecks << "\nWin pay: " << winPay << ":1\nBlackjack pay: " << blackjackPay << ":1\n" << std::endl;
	bool _continue = Helper::getYesNo("Continue with options? (YES/NO)");
	if (_continue) {
		// create empty vector for human players
		// add one to factors to include initial bet
		winPay++;
		blackjackPay++;

		std::vector<std::string> playerNames;
		std::pair<std::vector<HumanPlayer>, std::vector<ComputerPlayer>> pair = getPlayerVector(0, 2, playerNames, 0);
		Blackjack* blackjack = new Blackjack(pair.first, pair.second, nDecks, winPay, blackjackPay);
		blackjack->runSimulation(nTimes);
	}
	else {
		runBlackjackSimulation();
	}
}

std::pair<std::vector<HumanPlayer>, std::vector<ComputerPlayer>> Menu::getPlayerVector(int nHuman, int nAI, std::vector<std::string> playerNames, float startingMoney) {

	std::vector<ComputerPlayer::PlayerType> aiType = { ComputerPlayer::NPC_ALEX,ComputerPlayer::NPC_SIDNEY };
	std::vector<std::string> genericNames = { "Alex ", "Sidney " };
	int nTypes = (int)aiType.size();

	std::vector<HumanPlayer> humanVec(nHuman);
	std::vector<ComputerPlayer> npcVec(nAI);

	for (int i = 0; i < nHuman; i++) {
		humanVec[i] = HumanPlayer(startingMoney, playerNames[i]);
	}

	for (int i = 0; i < nAI; i++) {
		std::string str = genericNames[i%nTypes] + std::to_string(i / nTypes + 1);
		// creates computer player, cycling between different ai types in 'aiType'
		npcVec[i] = ComputerPlayer(startingMoney, str, aiType[i%nTypes]);
	}
	std::pair<std::vector<HumanPlayer>, std::vector<ComputerPlayer>> pair = { humanVec,npcVec };
	return pair;
}