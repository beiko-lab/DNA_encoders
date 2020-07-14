class BaseForEncoder:
    def __init__(self):
        self.binary = None

    def synonimous(self, class_name):

        nonpromoter = ["non-promoter", "Sigma--", "Sigma00", "no-promoter", "nopromoter", "nonpromoter", 0, "0"]
        if self.binary:
            promoter = ["promoter", "Sigma++",
                        "Sigma19", "sigma19", 1, "1",
                        "Sigma24", "sigma24", 2, "2",
                        "Sigma28", "sigma28", 3, "3",
                        "Sigma32", "sigma32", 4, "4",
                        "Sigma38", "sigma38", 5, "5",
                        "Sigma54", "sigma54", 6, "6",
                        "Sigma70", "sigma70", 7, "7"]
        else:
            sigma19 = ["Sigma19", "sigma19", 1, "1"]
            sigma24 = ["Sigma24", "sigma24", 2, "2"]
            sigma28 = ["Sigma28", "sigma28", 3, "3"]
            sigma32 = ["Sigma32", "sigma32", 4, "4"]
            sigma38 = ["Sigma38", "sigma38", 5, "5"]
            sigma54 = ["Sigma54", "sigma54", 6, "6"]
            sigma70 = ["Sigma70", "sigma70", 7, "7"]

        if class_name in nonpromoter:
            return nonpromoter[0]
        elif class_name in nonpromoter:
            return promoter[0]
        else:
            print("BaseForEncoder.synonimous is not implemented for your class named "+str(class_name))
            exit()