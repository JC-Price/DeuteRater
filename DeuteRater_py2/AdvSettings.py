"""
Copyright (c) 2016 Bradley Naylor, Michael Porter, J.C. Price, and Brigham Young University

All rights reserved.

Redistribution and use in source and binary forms,
with or without modification, are permitted provided
that the following conditions are met:

    * Redistributions of source code must retain the
      above copyright notice, this list of conditions
      and the following disclaimer.
    * Redistributions in binary form must reproduce
      the above copyright notice, this list of conditions
      and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors
      may be used to endorse or promote products derived
      from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""



from PyQt4 import QtCore, QtGui, uic
import File_Alterer as FA
import os
location =  os.path.dirname(os.path.abspath(__file__))
form_class = uic.loadUiType(os.path.join(location,"Advanced_Settings.ui"))[0]


class Change_Adv_Settings(QtGui.QDialog, form_class):
    def __init__(self, parent = None, user = None, default = None, user_s = None, s_filters = None):
        super(Change_Adv_Settings, self).__init__(parent)
        self.setWindowTitle("Advanced Settings Menu")
        self.setupUi(self)
        self.current = user
        self.user_settings = user_s
        self.static_filters = s_filters
        self.elem_comp_button.clicked.connect(self.elemental_change)
        self.mod_button.clicked.connect(self.modification_change)
        self.label_button.clicked.connect(self.labeling_change)
        self.ExitButton.clicked.connect(self.close)
        self.ApplyButton.clicked.connect(self.apply)
        self.RecalcButton.clicked.connect(self.recalc)
        if not self.current[0]:
            self.SlowOption.toggle()
    
    def elemental_change(self):
        self.changed_a_file = True
        self.elem_comp_menu = FA.alter_element_composition()
        self.elem_comp_menu.show()
    def modification_change(self):
        self.changed_a_file = True
        self.mod_menu = FA.alter_mods()
        self.mod_menu.show()
    def labeling_change(self):
        self.changed_a_file = True
        self.mod_menu = FA.amount_of_labels()
        self.mod_menu.show()
    def recalc(self):
        import MIDA_Calculator_vrs4 as mc
        message = mc.Recalculate(self.user_settings, self.static_filters)
        QtGui.QMessageBox.information(self, "Done", message)  

    def apply(self):
        if self.SlowOption.isChecked():
            self.current[0] = False
        if self.FastOption.isChecked():
            self.current[0] = True    
