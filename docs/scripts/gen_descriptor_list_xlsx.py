from openpyxl import Workbook
from openpyxl.styles import Font
from gen_default_descs import DescriptorInfo, get_all_descriptors


def main(out):
    wb = Workbook()
    ws = wb.active
    ws.title = "descriptors"

    ws.append(["index", "module", "name", "constructor", "dimention", "description"])
    for cell in ws["1:1"]:
        cell.font = Font(bold=True)

    i = 0
    for descs in get_all_descriptors():
        first = True
        for desc in descs:
            i += 1
            info = DescriptorInfo(desc)
            ws.append([
                i,
                info.module.split(".")[-1] if first else None,
                str(info.raw),
                "{}({})".format(info.constructor, ", ".join(info.parameters)),
                info.dimention,
                info.description,
            ])
            first = False

    for cells in ws.columns:
        L = max(len(str(cell.value or "")) for cell in cells)
        ws.column_dimensions[cells[0].column].width = L + 1

    wb.save(out)


if __name__ == "__main__":
    import sys
    main(sys.argv[1])
