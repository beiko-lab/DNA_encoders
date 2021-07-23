class BaseForEncoder:
    def __init__(self, binary=True):
        self.binary = binary

    def synonyms(self, class_name):

        nonpromoter = ["non-promoter", "Sigma--", "Sigma00", "no-promoter", "nopromoter", "nonpromoter", 0, "0"]
        if class_name in nonpromoter:
            return nonpromoter[0]

        if self.binary:
            promoter = ["promoter", "Sigma++",
                        "Sigma19", "sigma19", 1, "1",
                        "Sigma24", "sigma24", 2, "2",
                        "Sigma28", "sigma28", 3, "3",
                        "Sigma32", "sigma32", 4, "4",
                        "Sigma38", "sigma38", 5, "5",
                        "Sigma54", "sigma54", 6, "6",
                        "Sigma70", "sigma70", 7, "7"]
            if class_name in promoter:
                return promoter[0]
        else:
            if class_name in ["promoter", "Sigma++"]:
                return "promoter"
            if class_name in ["Sigma19", "sigma19", 1, "1"]:
                return "Sigma19"
            if class_name in ["Sigma24", "sigma24", 2, "2"]:
                return "Sigma24"
            if class_name in ["Sigma28", "sigma28", 3, "3"]:
                return "Sigma28"
            if class_name in ["Sigma32", "sigma32", 4, "4"]:
                return "Sigma32"
            if class_name in ["Sigma38", "sigma38", 5, "5"]:
                return "Sigma38"
            if class_name in ["Sigma54", "sigma54", 6, "6"]:
                return "Sigma54"
            if class_name in ["Sigma70", "sigma70", 7, "7"]:
                return "Sigma70"

        print(
            "BaseForEncoder.synonyms is not implemented for your class named "
            + str(class_name)
        )
        exit()