#include "sx.h"
#include "file.h"
#include "vid.h"
#include "world.h"

#include <stdlib.h>


int main(int argc, char *argv[])
{
  // global init
  sx_init(argc, argv);

  const char *mission = "c1m1.mis";
  if(argc > 1 && argv[1][0] != '-') mission = argv[1];

  const char *filename = mission;
  FILE *f = file_open(filename);
  strncpy(sx.mission.filename, filename, sizeof(sx.mission.filename));
  c3_mission_load(&sx.mission, f);
  fclose(f);

  sx_vid_start_mission();
  c3_mission_begin(&sx.mission);
  uint32_t last_event = SDL_GetTicks();
  uint32_t frames = 0;
  sx.time = last_event;
  uint32_t sim_time = last_event;
  const uint32_t delta_sim_time = 1000.0f/60.0f; // usual vsync
  int paused = 0;

  /*const char *vertical_logname = "/media/matze/T7 SSD/Bachelor/vertical_speed.txt";
  FILE *vertical = fopen(vertical_logname, "w");
  const char *horizontal_logname = "/media/matze/T7 SSD/Bachelor/horizontal_speed.txt";
  FILE *horizontal = fopen(horizontal_logname, "w");
  if (vertical == NULL || horizontal == NULL) 
  {
    printf("Error opening log file!\n");
  }
  else
  {
    printf("Logs successfully opened!\n");
    fprintf(vertical, "Vertikale Geschwindigkeit in m/s:Zeit in s\n");
    fprintf(horizontal, "Horizontale Geschwindigkeit in m/s:Zeit in s\n");
  }*/

  while(1)
  {
    sx_vid_render_frame_rect();

    if(sx.paused && !paused)
    {
      paused = 1;
    }
    if(paused && !sx.paused)
    {
      last_event = sim_time = SDL_GetTicks();
      paused = 0;
    }
    if(paused)
    {
      if(sx_vid_handle_input()) goto out;
      continue;
    }

    uint32_t end = SDL_GetTicks();
    while(sim_time < end)
    {
      if(sx_vid_handle_input()) goto out;
      // TODO: read network input and write movement controllers + rigid bodies
      // sx_net_recv();
      sx_world_think(delta_sim_time);
      sx_world_move(delta_sim_time);
      sim_time += delta_sim_time;
      // sx_net_send();
      // TODO: send out our movement controllers and rigid body locations with
      // current timestamp
    }
    // another second passed, run game mechanics only then
    frames++;
    if(end - last_event > 1000)
    {
      c3_mission_pump_events(&sx.mission);
      fprintf(stderr, "\r %.3g ms", (end - last_event)/(double)frames);
      frames = 0;
      last_event = end;
    }
    float horizontal_speed[3] = {0.0f, 0.0f, 0.0f};
    horizontal_speed[0] = sx.world.entity[sx.world.player_entity].body.v[0];
    horizontal_speed[2] = sx.world.entity[sx.world.player_entity].body.v[2];
    //fprintf(vertical, "%f:%f\n", sx.world.entity[sx.world.player_entity].body.v[1], sx.time / 1000.0f);
    //fprintf(horizontal, "%f:%f\n", length(horizontal_speed), sx.time / 1000.0f);
    // time is real time (lower gear animation etc)
    sx.time = end;
  }
out:
  /*fclose(vertical);
  fclose(horizontal);*/
  sx_vid_end_mission();
  c3_mission_end(&sx.mission);
  sx_cleanup(); // also cleans vid
  exit(0);
}
