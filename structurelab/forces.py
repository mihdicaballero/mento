class Forces:
    def __init__(self, My=0.0):
        self.My = My

    def get_moment(self):
        return self.My

    def set_moment(self, My):
        self.My = My