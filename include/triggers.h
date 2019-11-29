#pragma once
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

// general layout is one line per command, like:
// if
// condition
// then
// action

// TODO:
// cond
// sees attacked
// action
// setai
// route

typedef enum c3_condition_midi_t
{
  C3_COND_DIRE = 0,
  C3_COND_FLIGHT,
  C3_COND_LOSE_WIN,
  C3_COND_WIN,
  C3_COND_LOSE,
  C3_COND_COMBAT,
  C3_COND_PAD,
  C3_COND_SCARY,
  C3_COND_MIDI_SIZE,
}
c3_condition_midi_t;


static const char *c3_condition_midi_text[] =
{
  "dire.mid\0",
  "flight.mid\0",
  "lwin.mid\0",
  "win.mid\0",
  "lose.mid\0",
  "combat.mid\0",
  "pad.mid\0",
  "scary.mid\0",
};

typedef enum c3_condition_t
{
  C3_COND_NEARER = 0, // nearer, <distance>, <object/waypoint>
  C3_COND_FARTHER,    // farther, <distance>, <object>
  C3_COND_ALIVE,      // alive, <object>
  C3_COND_KILLED,     // killed, <object>, <num??>
  C3_COND_DESTROYED,  // killed, <object>
  C3_COND_BELOW,      // below, <altitude>
  C3_COND_ABOVE,      // above, <altitude>
  C3_COND_COUNTER,    // counter, <value>
  C3_COND_TIME,       // time, <value>
  C3_COND_WAYPOINT,   // waypoint, <object>, <waypoint id>
  C3_COND_DAMAGE,     // damage, <num>
  C3_COND_ENCOUNTER,  // encounter, <object>, <object>
  C3_COND_FACES,      // faces, <object/waypoint id>  // maybe looking at it?
  C3_COND_WEAPON,     // weapon, <num> // have this weapon selected
  C3_COND_CNT,
}
c3_condition_t;

static const char *c3_condition_text[] =
{
  "nearer",
  "farther",
  "alive",
  "killed",
  "destroyed",
  "below",
  "above",
  "counter",
  "time",
  "waypoint",
  "damage",
  "encounter",
  "faces",
  "weapon",
  "---",
};

typedef enum c3_operator_t
{
  C3_OP_AND = 0,
  C3_OP_ANDNOT,
  C3_OP_OR,
  C3_OP_CNT,
}
c3_operator_t;

static const char *c3_operator_text[] =
{
  "and",
  "andN",
  "or",
  "---",
};

typedef enum c3_action_t
{
  C3_ACT_TEXT = 0,     // text, <num>  // show text from inf file
  C3_ACT_FLASH,        // flash, <id>, <times>  // tutorial: flash HUD element
  C3_ACT_PLAY,         // play <wav filename>
  C3_ACT_VAPORIZE,     // vaporize, <object> ???
  C3_ACT_ELIMINATE,    // eliminate, <object> ???
  C3_ACT_WIN,          // mission won
  C3_ACT_CLEARCOUNTER, // restart from 0?
  C3_ACT_CNT,
}
c3_action_t;

static const char *c3_action_text[] =
{
  "text",
  "flash",
  "play",
  "vaporize",
  "eliminate",
  "win",
  "clearcounter",
  "---",
};

typedef struct c3_trigger_t
{
  c3_condition_t cond[32];
  uint32_t       arg0[32];
  uint32_t       arg1[32];
  c3_operator_t  oper[32];
  c3_action_t    action;
  uint32_t       aarg0;
  uint32_t       aarg1;
  uint32_t       enabled;
}
c3_trigger_t;

// forward declare
typedef struct c3_mission_t c3_mission_t;

int c3_triggers_parse(c3_mission_t *mis, FILE *f);

void c3_trigger_action(c3_mission_t *mis, c3_action_t a, uint32_t arg0, uint32_t arg1);

int c3_triggers_check_cond(c3_mission_t *mis, c3_condition_t c, uint32_t arg0, uint32_t arg1);

void c3_triggers_check_one(c3_mission_t *mis, c3_trigger_t *t);

void c3_triggers_check(c3_mission_t *mis);

char* c3_triggers_parse_music(char* filename, char letter, int gamestate, char fg);
