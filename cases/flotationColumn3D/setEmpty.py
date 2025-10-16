

# Открываем файл для чтения
with open('constant/polyMesh/boundary', 'r') as file:
    lines = file.readlines()

# Проходим по строкам и заменяем 'patch' на 'empty' в блоке frontAndBack
modified_lines = []
in_front_and_back_block = False

for line in lines:
    tmp = line
    if 'defaultFaces' in line:
        in_front_and_back_block = True

    if in_front_and_back_block:
        if 'type            empty;' in line:
            tmp = line.replace('empty', 'wall')
        if 'inGroups        1(empty);' in line:
            tmp = line.replace('empty', 'wall')
        
    if in_front_and_back_block and '}' in line:
        in_front_and_back_block = False

    modified_lines.append(tmp)

# Сохраняем изменения обратно в тот же файл
with open('constant/polyMesh/boundary', 'w') as file:
    file.writelines(modified_lines)

