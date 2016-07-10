#include "Blackjack.h"

template int Blackjack::processMoves<HumanPlayer>(GenericPlayer::Moves, HumanPlayer&, int);
template int Blackjack::processMoves<ComputerPlayer>(GenericPlayer::Moves, ComputerPlayer&, int);
template int Blackjack::processMoves<Dealer>(GenericPlayer::Moves, Dealer&, int);

template Blackjack::Result Blackjack::processResults<HumanPlayer>(HumanPlayer&);
template Blackjack::Result Blackjack::processResults<ComputerPlayer>(ComputerPlayer&);

template Blackjack::Result Blackjack::getResult<HumanPlayer>(HumanPlayer&);
template Blackjack::Result Blackjack::getResult<ComputerPlayer>(ComputerPlayer&);

Blackjack::Blackjack(std::vector<HumanPlayer> playerList, std::vector<ComputerPlayer> npcList, int nDecks, bool print, bool includeBets) {
	this->dealer = Dealer(nDecks);
	this->playerList = playerList;
	this->npcList = npcList;
	this->print = print;
	this->includeBets = includeBets;
}

Blackjack::Blackjack(std::vector<HumanPlayer> playerList, std::vector<ComputerPlayer> npcList, int nDecks, float winPayFactor, float blackjackPayFactor) {
	this->print = false;
	this->includeBets = false;
	this->dealer = Dealer(nDecks);
	this->playerList = playerList;
	this->npcList = npcList;
	this->winPayFactor = winPayFactor;
	this->blackjackPayFactor = blackjackPayFactor;
}

void Blackjack::run() {
	// runs a game of blackjack using the list of human and npc players
	std::cout << "======================Playing Blackjack======================" << std::endl;

	do {
		// keep playing hands until there are no user players 
		// if there are no human players (because they quit / run out of money or the user chose to have 0 human players), then 
		// ask every hand it should continue
		newHand();
		handNumber++;
	} while (playerList.size() > 0 || (Helper::getYesNo("New Hand?") && npcList.size()>0));
	std::cout << "Blackjack finished" << std::endl;
	this->~Blackjack();
}

void Blackjack::runSimulation(int nTimes) {
	// four results for each hand; blackjack, win, lose, push
	std::cout << "======================Blackjack simulation======================" << std::endl;
	std::vector<std::vector<int>> results;

	for (unsigned int i = 0; i < npcList.size(); i++) {
		// results is a 2d vector corresponding to the amount of times each player gets one of the four results;
		// int the inner vector, each index corresponds to each value in the enum 'Result' (ie blackjack win lose push)
		std::vector<int> _newVec = { 0,0,0,0 };
		results.push_back(_newVec);
	}
	for (int n = 0; n < nTimes; n++) {
		//running simulation nTimes
		std::vector<Result> handResultsVector = newHand();
		for (unsigned int i = 0; i < npcList.size(); i++) {
			results[i][handResultsVector[i]]++;
		}
	}
	for (unsigned int i = 0; i < npcList.size(); i++) {
		//printing results:
		//for profit, assume player bets $1 each hand; each win would give him back $2, a blackjack $2.50 etc
		float percentageProfit = (winPayFactor * results[i][WIN] + blackjackPayFactor*results[i][BLACKJACK] + results[i][PUSH] - nTimes) / ((float)nTimes) * 100;
		std::cout << "=======\n" << npcList[i].getName() << "\n=======" << std::endl;
		std::cout << "Win percentage (including blackjacks) : " << (results[i][BLACKJACK] + results[i][WIN]) / (float)nTimes * 100 << "%" << std::endl;
		std::cout << "Loss percentage : " << (results[i][LOSE]) / (float)nTimes * 100 << "%" << std::endl;
		std::cout << "Blackjack percentage : " << (results[i][BLACKJACK]) / (float)nTimes * 100 << "%" << std::endl;
		std::cout << "Percentage profit : " << percentageProfit << "%" << std::endl;
	}

}

std::vector<Blackjack::Result> Blackjack::newHand() {
	float minBet = 2;
	float maxBet = 20;
	if (includeBets) {
		makeBets(minBet, maxBet);
	}
	if (dealer.getDeckSize() < 40) {
		dealer.newDeck();
	}
	initialDeal();
	for (unsigned int i = 0; i < playerList.size(); i++) {
		i += processMoves(playerList[i].getMove(print), playerList[i], i);
	}
	for (unsigned int i = 0; i < npcList.size(); i++) {
		processMoves(npcList[i].getMove(print), npcList[i], 0);
	}
	if (print) {
		dealer.printHand();
	}
	processMoves(dealer.getMove(print), dealer, 0);

	return handResults();
}

void Blackjack::initialDeal() {
	// deals 2 cards to each player
	for (unsigned int i = 0; i < playerList.size(); i++) {
		playerList[i].clearHand();
		dealer.dealCards<HumanPlayer>(2, playerList[i], print);

	}

	for (unsigned int i = 0; i<npcList.size(); i++) {
		npcList[i].clearHand();
		dealer.dealCards<ComputerPlayer>(2, npcList[i], print);
	}

	dealer.clearHand();
	dealer.dealSelf(print);
}

void Blackjack::makeBets(float minBet, float maxBet) {
	if (print) {
		std::cout << "======================Hand " << handNumber << " bets======================" << std::endl;
	}

	for (unsigned int i = 0; i < playerList.size(); i++) {
		if (playerList[i].getMoney() >= minBet) {
			float max = Helper::getMin(playerList[i].getMoney(), maxBet);
			if (print) {
				std::cout << "Player " << playerList[i].getName() << " money: " << playerList[i].getMoney() << std::endl;
			}
			std::string str2 = "Enter bet for player " + playerList[i].getName() + " (max bet " + std::to_string(maxBet) + ", min bet " + std::to_string(minBet) + ")";
			float bet = Helper::getFloat(minBet, max, str2);
			playerList[i].setBet(bet);
			if (print) {
				std::cout << "Player " << playerList[i].getName() << " bets " << playerList[i].getBet() << ", current money: " << playerList[i].getMoney() << std::endl;
				std::cout << "------------------" << std::endl;
			}
		}
		else {
			playerList.erase(playerList.begin() + i);
			i--;
			if (print) {
				std::cout << "Player " << playerList[i].getName() << " cannot afford minimum bet and is out" << std::endl;
				std::cout << "------------------" << std::endl;
			}
		}
	}

	for (unsigned int i = 0; i < npcList.size(); i++) {
		if (npcList[i].getMoney() >= minBet) {
			npcList[i].setBet(minBet);
			if (print) {
				std::cout << "Player " << npcList[i].getName() << " bets " << npcList[i].getBet() << ", current money: " << npcList[i].getMoney() << std::endl;
				std::cout << "------------------" << std::endl;
			}
		}
		else {
			npcList.erase(npcList.begin() + i);
			i--;
			if (print) {
				std::cout << "Player " << npcList[i].getName() << " cannot afford minimum bet and is out" << std::endl;
				std::cout << "------------------" << std::endl;
			}
		}
	}

	if (print) {
		std::cout << "======================Hand " << handNumber << " deal======================" << std::endl;
	}
}

std::vector<Blackjack::Result> Blackjack::handResults() {
	// calls 'processResults' for each player; returns a vector for Sam-Alex test
	int n1 = playerList.size();
	int n2 = npcList.size();
	std::vector<Result> returnedVector;

	for (int i = 0; i < n1; i++) {
		returnedVector.push_back(processResults(playerList[i]));
	}

	for (int i = 0; i < n2; i++) {
		returnedVector.push_back(processResults(npcList[i]));
	}
	return returnedVector;
}

template <class T>
int Blackjack::processMoves(GenericPlayer::Moves moves, T &player, int index) {
	if (moves == GenericPlayer::HIT) {
		do {
			dealer.dealCards<T>(1, player, print);
		} while (player.getMove(print) == GenericPlayer::HIT);
	}
	if (moves == GenericPlayer::STICK) {
		return 0;
	}
	if (moves == GenericPlayer::DOUBLE) {
		player.doubleDown(print);
		dealer.dealCards(1, player, print);
	}
	if (moves == GenericPlayer::QUIT) {
		playerList.erase(playerList.begin() + index);
		return -1;
		//returns -1 if the user quits to account for playerList becoming shorter
	}
	if (print) {
		std::cout << "------------------" << std::endl;
	}

	return 0;
}

template <class T>
Blackjack::Result Blackjack::processResults(T &player) {
	Result result = getResult(player);
	if (print) {
		player.printHand();
	}
	if (result == LOSE) {
		// player gets no money returned. 
		// bet cleared for next round.
		player.setBet(0);
		if (print) {
			std::cout << "Player " << player.getName() << " loses" << std::endl;
		}
	}
	if (result == WIN) {
		// player gets original bet back plus 1:1 winnings
		player.addMoney(winPayFactor * player.getBet());
		if (print) {
			std::cout << "Player " << player.getName() << " beat the house and receives " << winPayFactor*player.getBet() << std::endl;
		}
		player.setBet(0);
	}
	if (result == PUSH) {
		player.addMoney(player.getBet());
		if (print) {
			std::cout << "Player " << player.getName() << " pushes and receives initial bet back" << std::endl;
		}
		player.setBet(0);
	}
	if (result == BLACKJACK) {
		player.addMoney(blackjackPayFactor*player.getBet());
		if (print) {
			std::cout << "Player " << player.getName() << " has a blackjack and receives " << blackjackPayFactor*player.getBet() << std::endl;
		}
		player.setBet(0);
	}
	return result;
}

template <class T>
Blackjack::Result Blackjack::getResult(T &player) {
	// returns different possible results given dealer and players hand
	int score = player.getScore();
	int dealerScore = dealer.getScore();
	bool blackjack = player.hasBlackjack();
	bool dealerBlackjack = dealer.hasBlackjack();
	if (score>21) {
		// if bust, automatic loss
		return LOSE;
	}
	if (dealerBlackjack && !blackjack) {
		// if dealer has blackjack and player doesn't, player loses
		return LOSE;
	}
	if (blackjack) {
		if (!dealerBlackjack) {
			// if the player has a blackjack and dealer doesn't, 3:2 payout
			return BLACKJACK;
		}
		else {
			// if both have blackjack, draw (push)
			return PUSH;
		}
	}
	if (score<21 && score == dealerScore) {
		// if score is same as dealers (and not bust), push
		// case of blackjacks covered earlier
		return PUSH;
	}
	if (score > dealerScore) {
		if (score <= 21) {
			return WIN;
		}
	}
	if (dealerScore > 21 && score <= 21) {
		// if dealer busts and player doesn;t
		return WIN;
	}
	return LOSE;
}
