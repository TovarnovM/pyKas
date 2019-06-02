import sqlite3
import random as rnd
from functional import seq
import datetime

class OptiSaver(object):
    """
    класс необходим для сохраниения/загрузки логов и непосредственно самих OptiPerson
    """
    default_db_path = "optisaver_tst.db"

    def __init__(self, conn):
        """
        Конструктор
        :param conn: ссл=ылка на sqlite3.connect(...)
        """
        self.conn = conn  # ссылка на sqlite3.connect(...)
        self.cursor = conn.cursor()  # ссылка на курсор
        self.create_main_table()  # создаем каркас в базе данных
        self.proj_name = None  # имя открытого проетка
        self.log_name = None  # имя таблицы для логов
        self.pop_name = None  # имя таблицы для хранения OptiPerson
        self.seed = None  # seed для рандомизатора

    @classmethod
    def open_project(cls, proj_name, db_path=None):
        """
        Классовый метод для открытия существующего проекта
        :param proj_name: имя проекта
        :param db_path: путь до базы данных (None _ для загрузки стандартного пути)
        :return: экземпляр OptiSaver

        :raise ValueError: если проект не найден
        """
        if not db_path:
            db_path = cls.default_db_path
        conn = sqlite3.connect(db_path)
        new_saver = cls(conn)
        if new_saver.open_proj(proj_name):
            return new_saver
        else:
            raise ValueError(f"{proj_name} не существует")

    @classmethod
    def create_project(cls, proj_name, db_path=None, seed=None):
        """
        Классовый метод для создания нового проекта
        :param proj_name: имя проекта
        :param db_path: путь до базы данных
        :param seed: seed для рандомизатора
        :return: экземпляр OptiSaver

        :raise ValueError: если проект уже существует
        """
        if not db_path:
            db_path = cls.default_db_path
        conn = sqlite3.connect(db_path)
        new_saver = cls(conn)
        if new_saver.has_proj(proj_name):
            raise ValueError(f"{proj_name} уже существует")
        new_saver.create_new_project(proj_name, seed)
        new_saver.open_proj(proj_name)
        return new_saver

    def has_proj(self, proj_name):
        """
        Проверяет есть ли проект в базе данных
        :param proj_name: имя проекта
        :return: True - если проект есть
        """
        lst = seq.sqlite3(self.conn, f"SELECT name FROM projects WHERE name='{proj_name}'").to_list()
        return len(lst) != 0

    def get_all_proj(self):
        """
        Получить список всех проектов
        :return: список всех проектов
        """
        lst = seq.sqlite3(self.conn, f"SELECT name FROM projects").map(lambda x: x[0]).to_list()
        return lst

    def open_proj(self, proj_name):
        """
        открывает (и записывает в поля объекта) существующий проект
        :param proj_name: имя проекта
        :return: True - если удалось открыть проект
        """
        if not self.has_proj(proj_name):
            return False
        self.proj_name, self.pop_name, self.log_name, self.seed = seq \
            .sqlite3(self.conn,
                     f"SELECT name,name_pop_table,name_log_table,rnd_seed FROM projects WHERE name='{proj_name}';") \
            .first()
        return True

    def close(self):
        """
        закрывает проект
        :return: none
        """
        self.conn.close()

    def create_new_project(self, proj_name, seed=None):
        """
        Создает новый проект
        :param proj_name: имя нового проекта
        :param seed: seed для рандомизатора
        :return: None
        """
        name_pop_table = "pop_" + proj_name
        name_log_table = "log_" + proj_name
        if not seed:
            seed = rnd.randint(0, 1e9)

        # добавляем строку в projects таблицу
        sql_command = f"""INSERT INTO projects (name, name_pop_table, name_log_table, date, rnd_seed)
            VALUES ("{proj_name}", "{name_pop_table}", "{name_log_table}", "{datetime.datetime.now()}", {seed});"""
        self.cursor.execute(sql_command)

        # создаем новую таблицу для логов
        sql_command = f"""
            CREATE TABLE IF NOT EXISTS {name_log_table} ( 
            id INTEGER PRIMARY KEY, 
            date DATE,
            generation INTEGER, 
            message TEXT);"""
        self.cursor.execute(sql_command)

        # создаем новую таблицу для сохранения популяций
        sql_command = f"""
            CREATE TABLE IF NOT EXISTS {name_pop_table} ( 
            id INTEGER PRIMARY KEY, 
            date DATE,
            generation INTEGER, 
            person TEXT);"""
        self.cursor.execute(sql_command)
        self.conn.commit()

    def delete_proj(self, proj_name):
        """
        Удаляет проект из базы данных
        :param proj_name: имя проекта для удаления
        :return: True - если успешно удалили проект
        """
        if not self.has_proj(proj_name):
            return False
        try:
            proj_id, name_pop, name_log = seq \
                .sqlite3(self.conn, f"SELECT id,name_pop_table,name_log_table FROM projects WHERE name='{proj_name}';") \
                .first()
            sql_command = f"DELETE FROM projects WHERE id={proj_id};"
            self.cursor.execute(sql_command)

            sql_command = f"DROP TABLE {name_pop};"
            self.cursor.execute(sql_command)

            sql_command = f"DROP TABLE {name_log};"
            self.cursor.execute(sql_command)

            self.conn.commit()
            return True
        except:
            return False

    def create_main_table(self):
        """
        Создает главную таблицу с проектами
        :return:
        """
        sql_command = """
            CREATE TABLE IF NOT EXISTS projects ( 
            id INTEGER PRIMARY KEY, 
            name TEXT UNIQUE, 
            name_pop_table TEXT, 
            name_log_table TEXT, 
            date DATE,
            rnd_seed INTEGER);"""
        self.cursor.execute(sql_command)
        self.conn.commit()

    def write_pops(self, gener_n, pop_list):
        """
        метод записи списка из string-ов в таблицу
        :param gener_n: номер поколения
        :param pop_list: список строк (сериализованных OptiPerson'ов)
        :return: None
        """
        date = datetime.datetime.now()
        for pop in pop_list:
            sql_command = f"""INSERT INTO {self.pop_name} 
                (date, generation, person)
                VALUES 
                ("{date}", {gener_n}, "{pop}");"""
            self.cursor.execute(sql_command)
        self.conn.commit()

    def delete_pops(self, gener_n):
        sql_command = f"""DELETE FROM {self.pop_name}
            WHERE generation={gener_n};
        """
        self.cursor.execute(sql_command)
        self.conn.commit()

    def write_log(self, gener_n, log_message):
        """
        Записать 1 запись лога
        :param gener_n: номер поколения
        :param log_message: лог-сообщение
        :return:
        """
        date = datetime.datetime.now()
        sql_command = f"""INSERT INTO {self.log_name} 
            (date, generation, message)
            VALUES 
            ("{date}", {gener_n}, "{log_message}");"""
        self.cursor.execute(sql_command)
        self.conn.commit()

    def load_generation(self, n_gener=None):
        """
        Загрузить список строк
        :param n_gener: номер поколения, из которого нужно загрузить поколение
        :return: список строк
        """
        if not n_gener:
            n_gener = self.get_last_generation()
        lst = seq.sqlite3(self.conn, f"""SELECT person 
            FROM {self.pop_name} 
            WHERE generation={n_gener}""") \
            .map(lambda x: x[0]) \
            .to_list()
        return lst

    def get_last_generation(self):
        """
        получить номер последнего сохраненного поколения
        :return: int поколения
        """
        return seq.sqlite3(self.conn,
                           f"""SELECT MAX(generation)
                FROM {self.pop_name} """) \
            .map(lambda x: x[0]) \
            .first()

    def load_logs(self, n_gener=None):
        """
        Загрузить логи
        :param n_gener: номер поколения
        :return: список строк
        """
        sql_command = f"""SELECT message 
            FROM {self.log_name} """
        if n_gener:
            sql_command += f"WHERE generation={n_gener}"
        return seq.sqlite3(self.conn, sql_command) \
            .map(lambda x: x[0]) \
            .to_list()


def main():
    pr_name = 'new_proj1'
    saver = OptiSaver.create_project(pr_name, seed=7)
    saver.write_pops(0, ["11", "11", "11", "11", "11"])
    saver.write_pops(2, ["22", "22", "22", "22", "22"])
    saver.write_pops(3, ["33", "33", "33", "33", "33"])
    saver.write_log(4, "ssss")
    print(saver.load_generation())
    print(saver.load_generation(1))
    print(saver.load_logs())
    print(saver.seed)
    saver.delete_proj(pr_name)
    saver.close()


if __name__ == '__main__':
    main()
