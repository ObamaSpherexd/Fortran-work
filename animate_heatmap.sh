 #!/bin/bash
 # Script to create animated GIF from heatmap frames
 
 echo "Creating animated GIF..."
convert -delay 20 -loop 0 anim_frame_000.png anim_frame_*.png heat_animation.gif
 echo "Animation saved: heat_animation.gif"
 echo "To view: display heat_animation.gif"
