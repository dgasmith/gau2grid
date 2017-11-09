"""
Contains a simple codegen helper
"""


class CodeGen(object):
    def __init__(self, indent="    "):
        self.indent_lvl = 0
        self.indent_tab = indent
        self.data = []

    def indent(self, lvl=1):
        """
        Indents the code one or more levels
        """

        self.indent_lvl += lvl

    def dedent(self, lvl=1):
        """
        Dedents the code one or more levels
        """

        self.indent_lvl -= lvl
        if self.indent_lvl < 0:
            raise ValueError("Indent level is negative!")

    def write(self, line):
        shift = self.indent_lvl * self.indent_tab
        self.data.append(shift + line)

    def blankline(self):
        self.data.append("")

    def repr(self, combine="\n"):
        return combine.join(self.data)
