# Planner - класс, необходимый для сохранения результатов серии вычислений, и выдачи еще не выполненных

import sqlitedict
import os

class Planner(object):
    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.plan_dict_cache=None
        # self.sqldict = sqlitedict.SqliteDict(
        #     filename=filename, 
        #     autocommit=kwargs.get('autocommit', True))

    def set_plan(self, plan):
        plan_dict = {f'{i}': point for i, point in enumerate(plan)}
        self.set_plan_dict(plan_dict)

    def set_plan_dict(self, plan_dict):
        with sqlitedict.SqliteDict(filename=self.filename) as sqldict:
            sqldict['plan_dict'] = plan_dict
            sqldict.commit()
        self.plan_dict_cache=None
    
    @property
    def plan_dict(self):
        if self.plan_dict_cache is None:
            with sqlitedict.SqliteDict(filename=self.filename) as sqldict:
                self.plan_dict_cache = sqldict.get('plan_dict', None)
        return self.plan_dict_cache

    @property
    def plan(self):
        plan_dict = self.plan_dict
        if plan_dict:
            return [val for val in plan_dict.values()]
        return None

    def get_calced_ids(self):
        """Возвращает множество id, которые уже посчитаны
        
        Returns:
            set(str) -- set посчитанных id
        """
        with sqlitedict.SqliteDict(filename=self.filename) as sqldict:
            return { id_ for id_ in sqldict.keys() if id_ != 'plan_dict' }


    def get_uncalced(self):
        """Возвращает список кортежей (id, X) точек, еще не посчитанных
        X - точка в плане, id типа str ее id) 
        id нужно для сохранения результата
        
        Returns:
            list((str, type(X))) -- -__-
        """
        plan_dict = self.plan_dict
        all_ids = set(plan_dict.keys())
        calced_ids = self.get_calced_ids()
        uncalced = all_ids.difference(calced_ids)
        return [(uid, plan_dict[uid]) for uid in uncalced]

    def get_calced(self, include_id=True):
        """Получить список кортежей (id, X, Y) (или (X,Y), если include_id==False)
        id - str
        X - точка в плане
        Y - результат расчета (сохраненный ранее)
        
        Keyword Arguments:
            include_id {bool} -- нужно ли добавлять в начало кортежа id (default: {True})
        
        Returns:
            list(tuple) -- -__-
        """
        calced_ids = self.get_calced_ids()
        plan_dict = self.plan_dict
        with sqlitedict.SqliteDict(filename=self.filename) as sqldict:
            if include_id:
                return [(cid, plan_dict[cid], sqldict[cid]) for cid in calced_ids]
            else:
                return [(plan_dict[cid], sqldict[cid]) for cid in calced_ids]

    def save(self, id, value):
        if id in self.plan_dict:
            with sqlitedict.SqliteDict(filename=self.filename) as sqldict:
                sqldict[id] = value
                sqldict.commit()
        else:
            raise AttributeError(f'Такого ключа для сохранения нет {id}')

    def clean(self):
        with sqlitedict.SqliteDict(filename=self.filename) as sqldict:
            for key in sqldict:
                if key != 'plan_dict':
                    sqldict.pop(key)
            sqldict.commit()
    
    def delete_file(self):
        try:
            os.remove(self.filename)
            print(f'File {self.filename} removed')
        except Exception as e:
            print(e)

    
if __name__ == "__main__":
    planner = Planner('test.db')
    planner.set_plan([[1,2,3], [4,5,6], 7, 8, 9])
    plan = planner.plan_dict
    planner.save(1, 'wewe')
    planner.get_uncalced()
    planner.delete_file()

