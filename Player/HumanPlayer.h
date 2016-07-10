#ifndef _H_HUMANPLAYER_H_
#define _H_HUMANPLAYER_H_

#include "GenericPlayer.h"
#include "Helper.h"
#include <string>

class HumanPlayer : public GenericPlayer {
public:
	HumanPlayer(float money, std::string name);
	HumanPlayer():GenericPlayer(0,"Player") {}
	GenericPlayer::Moves getMove(bool print);
	std::vector<GenericPlayer::Moves> getAvailableMoves();
	std::string getAvailableMovesStr(std::vector<Moves> &allowedMoves);

};

#endif