#ifndef _H_BLACKJACK_H_
#define _H_BLACKJACK_H_

#include <vector>
#include "Dealer.h"
#include "Card.h"
#include "HumanPlayer.h"
#include "ComputerPlayer.h"

class Blackjack {
	
public:
	Blackjack(std::vector<HumanPlayer> playerList, std::vector<ComputerPlayer> npcList, int nDecks, bool print, bool includeBets);
	Blackjack(std::vector<HumanPlayer> playerList, std::vector<ComputerPlayer> npcList, int nDecks, float winPayFactor, float blackjackPayFactor);

	// declaring winnings variables here; these include the original bet
	float winPayFactor = 2;
	float blackjackPayFactor = 2.5f;
	void run();
	void runSimulation(int nTimes);

	const enum Result {
		BLACKJACK,WIN, PUSH, LOSE
	};
protected:
	std::vector<HumanPlayer> playerList;
	std::vector<ComputerPlayer> npcList;
	int handNumber = 1;
	Dealer dealer;
	bool print;
	bool includeBets;

	std::vector<Result> newHand();
	void initialDeal();
	void makeBets(float minBet, float maxBet);
	std::vector<Result> handResults();
	
	template <class T>
	int processMoves(GenericPlayer::Moves moves, T &player, int index);
	
	template <class T>
	Result processResults(T &player);
	
	template <class T>
	Result getResult(T &player);
};

template int Blackjack::processMoves<HumanPlayer>(GenericPlayer::Moves, HumanPlayer&, int);
template int Blackjack::processMoves<ComputerPlayer>(GenericPlayer::Moves, ComputerPlayer&, int);
template int Blackjack::processMoves<Dealer>(GenericPlayer::Moves, Dealer&, int);

template Blackjack::Result Blackjack::processResults<HumanPlayer>(HumanPlayer&);
template Blackjack::Result Blackjack::processResults<ComputerPlayer>(ComputerPlayer&);

template Blackjack::Result Blackjack::getResult<HumanPlayer>(HumanPlayer&);
template Blackjack::Result Blackjack::getResult<ComputerPlayer>(ComputerPlayer&);

#endif
