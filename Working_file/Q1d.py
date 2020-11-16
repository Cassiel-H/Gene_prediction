# Jiachen Huo
# 260779330

import pandas as pd


def main():
    # read the file into dataframe and rename the header
    predicted = pd.read_csv('1c_out.gff3', delimiter="\t", names=["Seq_id", "type", "feature", "start", "end"])
    actualread = pd.read_csv('Vibrio_vulnificus.ASM74310v1.37.gff3', delimiter="\t", \
                             names=["Seq_id", "type", "feature", "start", "end", "score", "strand"])

    # in actual dataframe, restrict to only the "+" strand and CDS feature types
    i = actualread[actualread['feature'] == 'CDS']
    actual = i[i['strand'] == '+']

    # change the actual dataframe's start and end as type int
    actual['start'] = actual['start'].astype(int)
    actual['end'] = actual['end'].astype(int)

    to_pre, to_actual = compare(predicted, actual)

    for i in to_pre:
        to_pre[i] /= len(actual)
        to_pre[i] = "{:.2%}".format(to_pre[i])
    for i in to_actual:
        to_actual[i] /= len(predicted)
        to_actual[i] = "{:.2%}".format(to_actual[i])

    print("The accuracy as fraction of annotated genes on the positive strand to predicted genes are: " + "\n" + str(
        to_pre) + "\n")
    print("The accuracy as fraction of predicted genes to annotated genes on the positive strand are: " + "\n" + str(
        to_actual))


def compare(data1, data2):

    to_data1 = {'Both_match': 0, 'Only_start': 0, 'Only_end': 0, 'no_match': 0}

    to_data2 = {'Both_match': 0, 'Only_start': 0, 'Only_end': 0, 'no_match': 0}

    for i in data1.index:
        for j in data2.index:
            if data1['Seq_id'][i] == data2['Seq_id'][j]:

                # if both start and end matches
                if int(data1['start'][i]) == int(data2['start'][j]) and int(data1['end'][i]) == int(data2['end'][j]):
                    # print("sfsdfsdf")
                    to_data1['Both_match'] += 1
                    to_data2['Both_match'] += 1

                # if only start matched
                elif int(data1['start'][i]) == int(data2['start'][j]) and int(data1['end'][i]) != int(data2['end'][j]):
                    # print("好的")
                    to_data1['Only_start'] += 1
                    to_data2['Only_start'] += 1

                # if only end matched
                elif int(data1['start'][i]) != int(data2['start'][j]) and int(data1['end'][i]) == int(data2['end'][j]):
                    # print("???")
                    to_data1['Only_end'] += 1
                    to_data2['Only_end'] += 1

                # if nothing matched
                else:
                    # print("what?")
                    to_data1['no_match'] += 1
                    to_data2['no_match'] += 1
                    break

    return to_data1, to_data2


if __name__ == '__main__':
    main()
