from fastapi import APIRouter
import os
from .api_config_utils import RecordStubRoute

RECORD_STUBS = os.environ['RECORD_STUBS'] == '1' if 'RECORD_STUBS' in os.environ else False


def get_router():
    if RECORD_STUBS:
        return APIRouter(route_class=RecordStubRoute)
    else:
        return APIRouter()
