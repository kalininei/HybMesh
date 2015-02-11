#Plugins directory
PlugDir = .

.SECONDEXPANSION:
ResTarget = ComGridRes
PlugTarget = TransferListWidget
FormsTarget = HistoryWindow ComGridMain UniteGridsDlg AddUnfCircDlg MoveRotateDlg \
	      GridViewOpt


all: $(ResTarget) $(PlugTarget) $(FormsTarget)

#Resources
$(ResTarget): $$@.qrc
	pyrcc4 $< -o $@_rc.py

#Plugins
$(PlugTarget): $(PlugDir)/$$@.ui
	pyuic4 $< -o $(PlugDir)/ui_$@.py


#Application windows
$(FormsTarget): $$@.ui
	pyuic4 $< -o ui_$@.py
