 #!/bin/bash
 echo "Creating animated GIF..."
convert -delay 7 -loop 0 anim_frame_wave000.png anim_frame_wave*.png wave_animation.gif
 echo "Animation saved: wave_animation.gif"
