
class FileUtils:
    """ Common file handling methods used by the pipelines """

    @staticmethod
    def write_output(content, output_filename):
        """Write the given text content to a file"""
        try:
            with open(output_filename, 'w') as out:
                out.write(content)
        except IOError:
            print('Cannot open filename starting "{}"'.format(output_filename))
            raise

    @staticmethod
    def create_output_contents(final_dict):
        """Create tab-delimited output file table from the given dictionary"""
        final = sorted(final_dict.items(), key=lambda item: item[0], reverse=False)
        content = ''
        for n, item in enumerate(final):
            if n == len(final)-1:
                content += item[0] + '\n'
            else:
                content += item[0] + '\t'
        for n, item in enumerate(final):
            if n == len(final)-1:
                content += item[1] + '\n'
            else:
                content += item[1] + '\t'
        return content

    @staticmethod
    def create_header_line(headers):
        output = 'ID'
        for header in headers:
            output += '\t' + header
        output += '\n'

        return output

    @staticmethod
    def read_delimited_id_file_with_hdrs(delimiter, input_filename, expected_num_fields):
        """Parse an input file into a dictionary keyed on the first column"""
        """Each dictionary 'value' is itself a dictionary of fields for that id key"""
        headers = []
        contents = {}
        try:
            with open(input_filename, 'r') as fd:
                lines = fd.readlines()
                line_num = 0
                for line in lines:
                    line_num += 1
                    line = line.strip()
                    if line == '':
                        continue

                    fields = line.split(delimiter)
                    if len(fields) != expected_num_fields:
                        raise RuntimeError('Invalid file format for ' + input_filename)

                    if line_num == 1:
                        # read column headers
                        headers = fields[1::]
                    else:
                        field_dict = {}
                        field_idx = 1
                        for header in headers:
                            field_dict[header] = fields[field_idx]
                            field_idx += 1
                        contents[fields[0]] = field_dict

        except IOError:
            print('Cannot open filename starting "{}"'.format(input_filename))
            raise

        return headers, contents

