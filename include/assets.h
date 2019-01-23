#pragma once

#include "sound.h"
#include "c3object.h"

#include <stdint.h>
#include <assert.h>
#include <string.h>


typedef struct sx_assets_t
{
  sx_sound_t sound[512];
  uint32_t num_sounds;

  sx_object_t object[512];  // these are the unique objects
  uint32_t num_objects;     // entities in world struct make use of these
}
sx_assets_t;

// returns handle and dedupes
// please only pass lower case strings
uint32_t sx_assets_load_sound(
    sx_assets_t *a,
    const char *filename);

// returns handle and dedupes
// use lower case file names
// loads .ai files
uint32_t sx_assets_load_object(
    sx_assets_t *a,
    const char *filename);
// TODO: routine to nuke the whole thing and start from scratch
