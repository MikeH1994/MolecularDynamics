#ifndef _H_COMPUTERPLAYER_H_
#define _H_COMPUTERPLAYER_H_
#include "GenericPlayer.h"
#include <string>

class ComputerPlayer : public GenericPlayer {
protected:
	PlayerType character;
public:
	ComputerPlayer(float money, std::string name, PlayerType character);
	ComputerPlayer() {}
	Moves getMove(bool print);
	Moves moveAlex();
	Moves moveSid();
	Moves moveDirac();
};


#endif