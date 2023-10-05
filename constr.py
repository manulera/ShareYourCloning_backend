from pydantic import BaseModel, constr


class Person(BaseModel):
    name: constr(pattern=r'^[A-Z]+$')
    age: int
