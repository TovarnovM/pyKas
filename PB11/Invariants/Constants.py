class Constants(object):
    def __init__(self, gamma, b=0):
        """
        gamma - показатель адиабаты
        b - коволюм
        """
        self.b = b
        self.g = [1 for i in range(10)]
        self.g[0] = gamma
        self.g[1] = 0.5 * (gamma - 1) / gamma
        self.g[2] = 0.5 * (gamma + 1) / gamma
        self.g[3] = 2 * gamma / (gamma - 1)
        self.g[4] = 2 / (gamma - 1)
        self.g[5] = 2 / (gamma + 1)
        self.g[6] = (gamma - 1) / (gamma + 1)
        self.g[7] = 0.5 * (gamma - 1)
        self.g[8] = 1 / gamma
        self.g[9] = gamma - 1