def wrapstring(string, linemax, delim):
    strings = string.split(' ')
    lines = [strings[0]]
    for s in strings[1:]:
        if len(lines[-1]) != 0 and len(lines[-1])+len(s) > linemax:
            lines.append(s)
        else:
            lines[-1] += ' '+s
    return delim.join(lines)
    #import re
    #return re.sub("(.{1,%u})\s(.{1,%u})" % (linemax,linemax), '\g<1>'+delim+'\g<2>', string)

s = "aaagfksdfkjghdfksja afgsdfgkjsdfghkaaa gfjdkshgfkjds gfjdsg f"

print wrapstring(s, 50, ' \\\n')
