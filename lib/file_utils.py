
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

