from pathlib import Path

# Application path
application = str(Path('dist/SeqAnalyse.app').absolute())

# Output filename
filename = str(Path('dist/SeqAnalyse.dmg').absolute())

# Volume format
format = 'UDRW'

# Volume name
volume_name = 'SeqAnalyse'

# Volume size
size = '1g'

# Icon locations
icon_locations = {
    'SeqAnalyse.app': (140, 120),
    'Applications': (400, 120)
}

# Window position and size
window_rect = ((100, 100), (500, 300))

# Background
background = 'builtin-arrow'

# Files to include
files = [application]

# Symlinks to create
symlinks = {'Applications': '/Applications'}

# Icon size
icon_size = 64

# Text size
text_size = 12
