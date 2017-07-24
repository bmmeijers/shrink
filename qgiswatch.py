# Copy and paste below code inside the Python interactive console of QGIS.
# When another process writes to the file in attachHandler
# QGIS will refresh the map pane

def draw():
    print "clearing cache and redrawing" 
    qgis.utils.iface.mapCanvas().setCachingEnabled(False)
    qgis.utils.iface.mapCanvas().resetCachedContent()
    qgis.utils.iface.mapCanvas().refresh()

def attachHandler():
    from PyQt4.QtCore import QFileSystemWatcher
    watcher = QFileSystemWatcher()
    watcher.addPath('/home/martijn/signal')
    watcher.fileChanged.connect(draw)
    print "attached watch handler"

attachHandler()
draw()
