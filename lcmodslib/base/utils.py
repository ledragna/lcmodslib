def readlistnum(x, nmax=None):
    result = []
    for part in x.split(','):
        if '-' in part:
            a, b = part.split('-')
            if a == "":
                a = 0
            else:
                a = int(a)
            if b == "":
                if nmax:
                    b = int(nmax)
                else:
                    raise NotImplementedError()
            else:
                b = int(b)
            result.extend(range(a, b + 1))
        else:
            a = int(part)
            result.append(a)
    return result