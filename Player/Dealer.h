#ifndef _H_DEALER_H_
#define _H_DEALER_H_

#include "Deck.h"
#include "GenericPlayer.h"
#include "HumanPlayer.h"
#include "ComputerPlayer.h"
#include "Card.h"
#include <vector>
#include <iostream>
#include <string>

class Dealer :public GenericPlayer {

protected:
	Deck deck;
	float pot;
	int nDecks;
public:
	Dealer();
	Dealer(int nDecks);
	void dealSelf(bool print);
	Card drawNext();
	GenericPlayer::Moves getMove(bool print);
	int getDeckSize();
	void newDeck();
	//void dealCards(int nCards, HumanPlayer player, bool print);
	//void dealCards(int nCards, ComputerPlayer player, bool print);
	template <class T>
	void dealCards(int nCards, T &player, bool print);
};

template void Dealer::dealCards<HumanPlayer>(int,HumanPlayer&,bool);
template void Dealer::dealCards<ComputerPlayer>(int, ComputerPlayer&, bool);
template void Dealer::dealCards<Dealer>(int,Dealer&,bool);

#endif
